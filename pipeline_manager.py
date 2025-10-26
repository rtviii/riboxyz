# Create a new file: ribctl/lib/npet/pipeline_manager.py

import os
import sys
sys.path.append("/home/rtviii/dev/riboxyz")
import json
import numpy as np
import pyvista as pv
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Callable, Optional

import open3d as o3d
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.landmarks.ptc_via_trna import PTC_location

# Import your core computation functions
# (This assumes kdtree_approach.py and alphalib.py are accessible)
from ribctl.lib.npet.kdtree_approach import (
    landmark_constriction_site,
    ribosome_entities,
    filter_residues_parallel,
    transform_points_to_C0,
    create_point_cloud_mask,
    transform_points_from_C0,
    clip_pcd_via_ashape,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    apply_poisson_reconstruction
)
from ribctl.lib.npet.alphalib import (
    cif_to_point_cloud,
    quick_surface_points,
    fast_normal_estimation,
    validate_mesh_pyvista
)
from ribctl import RIBETL_DATA

# --- ADD THIS DEBUGGING ---
print(f"DEBUG [PipelineManager]: Imported RIBETL_DATA = '{RIBETL_DATA}'")
# --- END DEBUGGING ---

if not RIBETL_DATA:
    raise EnvironmentError("RIBETL_DATA environment variable is not set or not imported correctly.")

RIBETL_DATA_PATH = Path(RIBETL_DATA)

# --- ADD THIS DEBUGGING ---
print(f"DEBUG [PipelineManager]: Resolved RIBETL_DATA_PATH = '{RIBETL_DATA_PATH.resolve()}'")
# --- END DEBUGGING ---

@dataclass
class StageDefinition:
    """A declarative definition of a pipeline stage."""
    name: str
    description: str
    run_method_name: str # Name of the method in PipelineManager to run
    default_params: Dict[str, Any] = field(default_factory=dict)
    
    # Function to get all artifacts this stage *might* produce
    # Takes (assets_dir, artifacts_dir) -> returns {artifact_name: path}
    get_artifact_paths: Callable[[Path, Path], Dict[str, Path]] = field(default_factory=lambda: lambda a, b: {})
    visualization_hint: str = 'json' # 'mesh', 'point_cloud', 'json'

class PipelineManager:
    """
    Orchestrates the NPET pipeline logic, decoupled from any UI.
    Manages state via the file system ("Artifact-First").
    """
    
    def __init__(self, rcsb_id: str):
        self.rcsb_id      = rcsb_id.upper()
        self.assets_dir   = RIBETL_DATA_PATH / self.rcsb_id
        self.artifacts_dir = self.assets_dir / "artifacts"
        self.artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        self.stages    = self._define_stages()
        self.stage_map = {s.name: s for s in self.stages}
        
        # This context holds *in-memory* data passed between stages *during a run*
        self.context: Dict[str, Any] = {}

    def get_stage_artifacts(self, stage_name: str) -> Dict[str, Path]:
        """Get the defined artifact paths for a given stage."""
        stage = self.stage_map[stage_name]
        return stage.get_artifact_paths(self.assets_dir, self.artifacts_dir)

    def check_stage_status(self, stage_name: str) -> bool:
        """Check if a stage is 'complete' by verifying its artifacts exist."""
        paths = self.get_stage_artifacts(stage_name).values()
        if not paths:
            return True # Stages without artifacts (like Setup) are always "complete"
        return all(p.exists() for p in paths)

# In ribctl/lib/npet/pipeline_manager.py

    def run_stage(self, stage_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run a single stage by name.
        This is the function our QThread will call.
        """
        stage_def = self.stage_map[stage_name]
        
        # 1. Clear old artifacts for this stage
        # *** FIX: Do NOT delete artifacts for the 'Setup' stage ***
        if stage_name != "Setup":
            for path in self.get_stage_artifacts(stage_name).values():
                if path.exists():
                    print(f"Clearing old artifact: {path.name}")
                    os.remove(path)
        # *** END FIX ***
                
        # 2. Ensure dependencies (context) are met
        self._load_context_for_stage(stage_def)
        
        # 3. Run the actual compute method
        method = getattr(self, stage_def.run_method_name)
        stage_params = {**stage_def.default_params, **params}
        
        print(f"Running stage '{stage_name}' with params: {stage_params}")
        result = method(stage_params)
        
        # 4. Update the in-memory context
        if result:
            self.context.update(result)
        
        print(f"Stage '{stage_name}' complete.")
        return result

    def _load_context_for_stage(self, stage: StageDefinition):
        """
        Super simple dependency loading.
        If a context key is missing, try to load it from a predecessor's artifact.
        """
        # --- Dependencies for "PTC Identification" ---
        if "cifpath" not in self.context:
            setup_artifacts = self.get_stage_artifacts("Setup")
            cif_path = setup_artifacts.get("cif")
            if cif_path and cif_path.exists():
                self.context["cifpath"] = cif_path
                self.context["ro"]      = RibosomeOps(self.rcsb_id) # Re-create RibosomeOps
            else:
                # If Setup hasn't even run, we can't do anything.
                # In a real run, Setup is always first.
                if stage.name != "Setup":
                     raise FileNotFoundError(f"CIF file not found. Run 'Setup' first.")
            
        # --- Dependencies for "Entity Filtering" (and many others) ---
        if "ptc_pt" not in self.context and stage.name not in ["Setup", "PTC Identification"]:
            ptc_artifact = self.get_stage_artifacts("PTC Identification").get("ptc_json")
            if ptc_artifact and ptc_artifact.exists():
                with open(ptc_artifact, 'r') as f:
                    self.context["ptc_pt"] = np.array(json.load(f)['location'])
            else:
                 raise FileNotFoundError(f"PTC JSON not found. Run 'PTC Identification' first.")

        if "constriction_pt" not in self.context and stage.name not in ["Setup", "PTC Identification", "Constriction Site"]:
            # This one is tricky, it's not a file. It's fetched.
            # Let's just fetch it if missing.
            try:
                self.context["constriction_pt"] = landmark_constriction_site(self.rcsb_id)
            except Exception as e:
                raise RuntimeError(f"Could not load constriction point: {e}")
        
        # Add more dependency loads as needed (e.g., loading filtered_points.npy)

    # ---------------------------------------------------------------------
    # STAGE DEFINITIONS
    # ---------------------------------------------------------------------
    
    def _define_stages(self) -> List[StageDefinition]:
        """
        This is the new "single source of truth" for the pipeline.
        We define all stages, their params, their artifacts, and their logic.
        """
        return [
            StageDefinition(
                name="Setup",
                description="Initialize pipeline and load structure data",
                run_method_name="_run_setup",
                get_artifact_paths=lambda assets, artifacts: {
                    "cif": assets / f"{self.rcsb_id}.cif"
                },
                visualization_hint='json'
            ),
            StageDefinition(
                name="PTC Identification",
                description="Identify Peptidyl Transferase Center location",
                run_method_name="_run_ptc_identification",
                get_artifact_paths=lambda assets, artifacts: {
                    # This path is from your 7K00 tree output
                    "ptc_json": assets / f"{self.rcsb_id}_PTC.json"
                },
                visualization_hint='json'
            ),
            StageDefinition(
                name="Constriction Site",
                description="Identify ribosome exit tunnel constriction site (fetches from API)",
                run_method_name="_run_constriction_identification",
                # This stage *produces* no file, just context.
                visualization_hint='json'
            ),
            StageDefinition(
                name="Alpha Shape (Exterior)",
                description="Generate exterior ribosome surface mesh",
                run_method_name="_run_alpha_shape",
                default_params={
                    "d3d_alpha": 200, "d3d_tol": 10, "d3d_offset": 3,
                    "kdtree_radius": 40, "max_nn": 60, "tangent_planes_k": 20,
                    "PR_depth": 6, "PR_ptweight": 4
                },
                get_artifact_paths=lambda assets, artifacts: {
                    # From your 7K00 tree output
                    "ashape_mesh": assets / f"{self.rcsb_id}_ALPHA_SHAPE.ply",
                    "ashape_mesh_ascii": assets / f"{self.rcsb_id}_ALPHA_SHAPE_ascii.ply",
                    "structure_ptcloud": artifacts / f"{self.rcsb_id}_structure_ptcloud.npy",
                    "surface_points": artifacts / f"{self.rcsb_id}_alpha_surface_points.npy",
                    "normals": artifacts / f"{self.rcsb_id}_alpha_normals.ply"
                },
                visualization_hint='mesh'
            ),
            StageDefinition(
                name="Entity Filtering",
                description="Filter atoms within tunnel cylinder region",
                run_method_name="_run_entity_filtering",
                default_params={"radius": 35, "height": 120},
                get_artifact_paths=lambda assets, artifacts: {
                    "filtered_points": artifacts / f"{self.rcsb_id}_filtered_points.npy"
                },
                visualization_hint='point_cloud'
            ),
            StageDefinition(
                name="Point Cloud Processing",
                description="Transform points, create mask, get interior",
                run_method_name="_run_point_cloud_processing",
                default_params={"voxel_size": 1, "atom_size": 2},
                get_artifact_paths=lambda assets, artifacts: {
                    "transformed_points": artifacts / f"{self.rcsb_id}_transformed_points.npy",
                    "interior_points": artifacts / f"{self.rcsb_id}_interior_points.npy",
                    "back_projected_points": artifacts / f"{self.rcsb_id}_back_projected.npy",
                },
                visualization_hint='point_cloud'
            ),
            StageDefinition(
                name="Clustering",
                description="Cluster tunnel points using DBSCAN",
                run_method_name="_run_clustering",
                default_params={"epsilon": 5.5, "min_samples": 600},
                get_artifact_paths=lambda assets, artifacts: {
                    "largest_cluster": artifacts / f"{self.rcsb_id}_largest_cluster.npy",
                    # Add paths for other saved clusters if needed
                },
                visualization_hint='point_cloud'
            ),
            StageDefinition(
                name="Refinement",
                description="Refine clusters with secondary DBSCAN",
                run_method_name="_run_refinement",
                default_params={"epsilon": 3.5, "min_samples": 175},
                get_artifact_paths=lambda assets, artifacts: {
                    "refined_cluster": artifacts / f"{self.rcsb_id}_refined_cluster.npy"
                },
                visualization_hint='point_cloud'
            ),
            StageDefinition(
                name="Surface Extraction",
                description="Extract surface points from refined cluster",
                run_method_name="_run_surface_extraction",
                default_params={"alpha": 2, "tolerance": 1, "offset": 2},
                get_artifact_paths=lambda assets, artifacts: {
                    "surface_points": artifacts / f"{self.rcsb_id}_surface_points.npy"
                },
                visualization_hint='point_cloud'
            ),
            StageDefinition(
                name="Normal Estimation",
                description="Estimate normals for surface reconstruction",
                run_method_name="_run_normal_estimation",
                default_params={"kdtree_radius": 10, "kdtree_max_nn": 15, "correction_tangent_planes_n": 10},
                get_artifact_paths=lambda assets, artifacts: {
                    "normal_estimated_pcd": artifacts / f"{self.rcsb_id}_normal_estimated_pcd.ply"
                },
                visualization_hint='mesh' # .ply can be point cloud, but let's show normals
            ),
            StageDefinition(
                name="Mesh Reconstruction",
                description="Reconstruct tunnel mesh using Poisson",
                run_method_name="_run_mesh_reconstruction",
                default_params={"depth": 6, "ptweight": 3},
                get_artifact_paths=lambda assets, artifacts: {
                    "npet_mesh": assets / f"{self.rcsb_id}_NPET_MESH.ply",
                    "npet_mesh_ascii": assets / f"{self.rcsb_id}_NPET_MESH_ascii.ply"
                },
                visualization_hint='mesh'
            ),
            # Add ValidationStage...
        ]

    # ---------------------------------------------------------------------
    # STAGE IMPLEMENTATIONS (Ported from your ..._stage.py files)
    # ---------------------------------------------------------------------

# In ribctl/lib/npet/pipeline_manager.py

# In ribctl/lib/npet/pipeline_manager.py

    def _run_setup(self, params: Dict) -> Dict[str, Any]:
        """Ported from SetupStage"""
        
        cif_path = self.get_stage_artifacts("Setup").get("cif")

        if not cif_path:
             raise FileNotFoundError(f"CIF path object is None! Check StageDefinition.")

        print(f"DEBUG [_run_setup]: Attempting to find file at: '{cif_path}'")
        print(f"DEBUG [_run_setup]: Absolute resolved path: '{cif_path.resolve()}'")

        # --- NEW DEBUG ---
        # Use glob, which is a more direct way to ask the OS
        import glob
        glob_path = str(cif_path.parent / "*.cif")
        print(f"DEBUG [_run_setup]: Globbing for: '{glob_path}'")
        all_cifs = glob.glob(glob_path)
        print(f"DEBUG [_run_setup]: Glob found these CIF files: {all_cifs}")
        # --- END NEW DEBUG ---

        if not cif_path.exists():
            parent_dir = cif_path.parent
            print(f"DEBUG [_run_setup]: Parent dir '{parent_dir}' exists: {parent_dir.exists()}")
            print(f"DEBUG [_run_setup]: Contents of {parent_dir}: {os.listdir(parent_dir)}")
            raise FileNotFoundError(f"CIF file not found at {cif_path}")
        
        ro = RibosomeOps(self.rcsb_id)
        
        return {
            "cifpath": cif_path,
            "ro": ro,
        }

    def _run_ptc_identification(self, params: Dict) -> Dict[str, Any]:
        """Ported from PTCIdentificationStage"""
        ptc_info = PTC_location(self.rcsb_id)
        ptc_pt = np.array(ptc_info.location)
        
        ptc_json_path = self.get_stage_artifacts("PTC Identification")["ptc_json"]
        with open(ptc_json_path, 'w') as f:
            f.write(ptc_info.model_dump_json(indent=2))
            
        return {
            "ptc_info": ptc_info,
            "ptc_pt": ptc_pt
        }

    def _run_constriction_identification(self, params: Dict) -> Dict[str, Any]:
        """Ported from ConstrictionIdentificationStage"""
        try:
            constriction_pt = landmark_constriction_site(self.rcsb_id)
            return {"constriction_pt": constriction_pt}
        except Exception as e:
            raise RuntimeError(f"Failed to fetch constriction site: {e}")

    def _run_alpha_shape(self, params: Dict) -> Dict[str, Any]:
        """Ported from alphalib.py"""
        cif_path = self.context["cifpath"]
        ro       = self.context["ro"]
        
        ptcloud_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["structure_ptcloud"]
        surface_pts_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["surface_points"]
        normals_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["normals"]
        mesh_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["ashape_mesh"]

        print("Extracting point cloud from CIF...")
        chains = ro.first_assembly_auth_asym_ids()
        ptcloud = cif_to_point_cloud(cif_path, chains, do_atoms=True)
        np.save(ptcloud_path, ptcloud)
        
        print("Calculating surface points (Delaunay 3D)...")
        surface_pts = quick_surface_points(
            ptcloud, params["d3d_alpha"], params["d3d_tol"], params["d3d_offset"]
        )
        np.save(surface_pts_path, surface_pts)

        print("Estimating normals...")
        normal_estimated_pcd = fast_normal_estimation(
            surface_pts, params["kdtree_radius"], params["max_nn"], params["tangent_planes_k"]
        )
        o3d.io.write_point_cloud(str(normals_path), normal_estimated_pcd)
        
        print("Applying Poisson Reconstruction...")
        apply_poisson_reconstruction(
            str(normals_path),
            mesh_path,
            recon_depth=params["PR_depth"],
            recon_pt_weight=params["PR_ptweight"],
        )
        
        # Post-processing: keep largest component
        mesh = pv.read(mesh_path)
        labeled = mesh.connectivity(largest=True)
        labeled.save(mesh_path)
        
        watertight = validate_mesh_pyvista(labeled)
        print(f"Alpha shape mesh generated. Watertight: {watertight}")
        
        return {"ashape_mesh": mesh_path, "watertight": watertight}

    def _run_entity_filtering(self, params: Dict) -> Dict[str, Any]:
        """Ported from EntityFilteringStage"""
        cifpath = self.context["cifpath"]
        ptc_pt  = self.context["ptc_pt"]
        constriction_pt = self.context["constriction_pt"]
        
        # Logic for tunnel_debris (hardcoded)
        tunnel_debris = {"3J7Z": ["a", "7"], "5GAK": ["z"], "5NWY": ["s"], "7A5G": ["Y2"], "9F1D": ["BK"]}
        skip_chains = tunnel_debris.get(self.rcsb_id, [])
        
        residues = ribosome_entities(self.rcsb_id, cifpath, "R", skip_chains)
        
        filtered_residues = filter_residues_parallel(
            residues, ptc_pt, constriction_pt, params["radius"], params["height"]
        )
        
        filtered_points = np.array([
            atom.get_coord()
            for residue in filtered_residues
            for atom in residue.child_list
        ])
        
        filtered_points_path = self.get_stage_artifacts("Entity Filtering")["filtered_points"]
        np.save(filtered_points_path, filtered_points)
        
        return {
            "filtered_residues": filtered_residues,
            "filtered_points": filtered_points
        }

    def _run_point_cloud_processing(self, params: Dict) -> Dict[str, Any]:
        """Ported from PointCloudProcessingStage"""
        filtered_points = self.context["filtered_points"]
        ptc_pt = self.context["ptc_pt"]
        constriction_pt = self.context["constriction_pt"]
        ashape_mesh_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["ashape_mesh"]

        # 1. Transform points to cylinder coordinate system
        transformed_points = transform_points_to_C0(filtered_points, ptc_pt, constriction_pt)
        
        # 2. Create voxel mask
        mask, (x, y, z) = create_point_cloud_mask(
            transformed_points,
            radius=self.stage_map["Entity Filtering"].default_params["radius"],
            height=self.stage_map["Entity Filtering"].default_params["height"],
            voxel_size=params["voxel_size"],
            radius_around_point=params["atom_size"],
        )
        
        # 3. Get "empty" voxel centers
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        empty_voxels_C0 = np.column_stack((X[~mask], Y[~mask], Z[~mask]))
        
        # 4. Transform "empty" points back to world coordinates
        empty_in_world_coords = transform_points_from_C0(empty_voxels_C0, ptc_pt, constriction_pt)
        
        # 5. Clip empty points by the exterior alpha shape
        ashape_mesh = pv.read(ashape_mesh_path)
        interior_points, _ = clip_pcd_via_ashape(empty_in_world_coords, ashape_mesh)

        # Save artifacts
        np.save(self.get_stage_artifacts("Point Cloud Processing")["transformed_points"], transformed_points)
        np.save(self.get_stage_artifacts("Point Cloud Processing")["interior_points"], interior_points)
        np.save(self.get_stage_artifacts("Point Cloud Processing")["back_projected_points"], empty_in_world_coords)

        return {"interior_points": interior_points}

    def _run_clustering(self, params: Dict) -> Dict[str, Any]:
        """Ported from ClusteringStage"""
        interior_points = self.context["interior_points"]
        
        db, clusters_container = DBSCAN_capture(
            interior_points, params["epsilon"], params["min_samples"]
        )
        
        largest_cluster_pts, cluster_id = DBSCAN_pick_largest_cluster(clusters_container)
        
        print(f"Picked largest cluster #{cluster_id} with {len(largest_cluster_pts)} points.")
        
        # Save artifacts
        np.save(self.get_stage_artifacts("Clustering")["largest_cluster"], largest_cluster_pts)
        # We could save all clusters here if needed
        
        return {"largest_cluster": largest_cluster_pts}

    def _run_refinement(self, params: Dict) -> Dict[str, Any]:
        """Ported from RefinementStage"""
        largest_cluster = self.context["largest_cluster"]
        
        db, clusters_container = DBSCAN_capture(
            largest_cluster, params["epsilon"], params["min_samples"]
        )
        
        refined_cluster_pts, cluster_id = DBSCAN_pick_largest_cluster(clusters_container)
        
        print(f"Picked largest refined cluster #{cluster_id} with {len(refined_cluster_pts)} points.")

        np.save(self.get_stage_artifacts("Refinement")["refined_cluster"], refined_cluster_pts)
        
        return {"refined_cluster": refined_cluster_pts}

    def _run_surface_extraction(self, params: Dict) -> Dict[str, Any]:
        """Ported from SurfaceExtractionStage"""
        refined_cluster = self.context["refined_cluster"]
        
        surface_points = quick_surface_points(
            refined_cluster, params["alpha"], params["tolerance"], params["offset"]
        )
        
        np.save(self.get_stage_artifacts("Surface Extraction")["surface_points"], surface_points)
        
        return {"surface_points": surface_points}

    def _run_normal_estimation(self, params: Dict) -> Dict[str, Any]:
        """Ported from NormalEstimationStage"""
        surface_points = self.context["surface_points"]
        
        normal_estimated_pcd = fast_normal_estimation(
            surface_points,
            params["kdtree_radius"],
            params["kdtree_max_nn"],
            params["correction_tangent_planes_n"],
        )
        
        pcd_path = self.get_stage_artifacts("Normal Estimation")["normal_estimated_pcd"]
        o3d.io.write_point_cloud(str(pcd_path), normal_estimated_pcd)
        
        return {"normal_estimated_pcd_path": pcd_path}

    def _run_mesh_reconstruction(self, params: Dict) -> Dict[str, Any]:
        """Ported from MeshReconstructionStage"""
        pcd_path = self.context["normal_estimated_pcd_path"]
        mesh_path = self.get_stage_artifacts("Mesh Reconstruction")["npet_mesh"]
        
        apply_poisson_reconstruction(
            str(pcd_path),
            mesh_path,
            recon_depth=params["depth"],
            recon_pt_weight=params["ptweight"],
        )
        
        # Post-processing: keep largest component
        mesh = pv.read(mesh_path)
        labeled = mesh.connectivity(largest=True)
        labeled.save(mesh_path)
        
        watertight = validate_mesh_pyvista(labeled)
        print(f"NPET mesh generated. Watertight: {watertight}")
        
        return {"npet_mesh": mesh_path, "npet_watertight": watertight}