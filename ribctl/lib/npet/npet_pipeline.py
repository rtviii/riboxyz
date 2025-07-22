from pathlib import Path
from typing import Dict, List, Optional, Any

import numpy as np

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.kdtree_approach import (
    filter_residues_parallel,
    landmark_constriction_site,
    landmark_ptc,
    ribosome_entities,
)
from ribctl.lib.npet.pipeline.setup_stage import SetupStage
from ribctl.lib.npet.pipeline.ptc_identification_stage import PTCIdentificationStage
from ribctl.lib.npet.pipeline.constriction_identification_stage import (
    ConstrictionIdentificationStage,
)
from ribctl.lib.npet.pipeline.exterior_mesh_stage import AlphaShapeStage
from ribctl.lib.npet.pipeline.landmark_identification_stage import (
    LandmarkIdentificationStage,
)
from ribctl.lib.npet.pipeline.entity_filtering_stage import EntityFilteringStage
from ribctl.lib.npet.pipeline.point_cloud_processing_stage import (
    PointCloudProcessingStage,
)
from ribctl.lib.npet.pipeline.clustering_stage import ClusteringStage, generate_clusters_pdb
from ribctl.lib.npet.pipeline.refinement_stage import RefinementStage
from ribctl.lib.npet.pipeline.surface_extraction_stage import SurfaceExtractionStage
from ribctl.lib.npet.pipeline.normal_estimation_stage import NormalEstimationStage
from ribctl.lib.npet.pipeline.mesh_reconstruction_stage import MeshReconstructionStage
from ribctl.lib.npet.pipeline.validation_stage import ValidationStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import (
    NPETProcessingTracker,
)

def generate_final_cluster_pdb(
    refined_cluster: np.ndarray,
    output_path: str,
    rcsb_id: str
) -> str:
    """
    Generate a PDB file for the final refined cluster used for mesh generation.
    
    Args:
        refined_cluster: The final cluster points array
        output_path: Where to write the PDB file
        rcsb_id: PDB ID for header
        
    Returns:
        str: Path to the generated PDB file
    """
    
    with open(output_path, 'w') as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    FINAL TUNNEL CLUSTER                    {rcsb_id}\n")
        pdb_file.write(f"TITLE     FINAL REFINED CLUSTER FOR {rcsb_id}\n")
        pdb_file.write("REMARK    Final cluster used for mesh generation\n")
        
        atom_num = 1
        connect_records = []
        
        print(f"Writing final cluster with {len(refined_cluster)} points")
        
        # Write each point as a CA atom in chain A, residue 1
        for point_idx, (x, y, z) in enumerate(refined_cluster):
            pdb_file.write(
                f"ATOM  {atom_num:5d}  CA  ALA A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
            )
            
            # Create bonds between consecutive atoms
            if point_idx > 0:
                connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")
            
            atom_num += 1
        
        # Write chain terminator
        pdb_file.write(f"TER   {atom_num:5d}      ALA A   1\n")
        
        # Write all CONECT records
        for connect in connect_records:
            pdb_file.write(connect)
            
        # Write end record
        pdb_file.write("END\n")
    
    print(f"Generated final cluster PDB: {output_path}")
    return output_path

def generate_cylinder_residues_json(
    rcsb_id: str,
    cif_path: str,
    ptc_coords: np.ndarray,
    constriction_coords: np.ndarray,
    output_path: str,
    radius: float = 35.0,
    skip_nascent_chain: List[str] = [],
) -> dict:
    """
    Generate cylinder residues JSON for molstar visualization.

    Args:
        rcsb_id: PDB ID
        cif_path: Path to mmCIF file
        ptc_coords: PTC coordinates from pipeline context
        constriction_coords: Constriction site coordinates from pipeline context
        output_path: Where to write the JSON file
        radius: Cylinder radius in Angstroms
        skip_nascent_chain: Chain IDs to skip

    Returns:
        dict: The cylinder residues data structure
    """
    import json

    # Define cylinder between PTC and constriction site
    base_point = ptc_coords
    axis_point = constriction_coords

    # Calculate height as distance between landmarks plus some buffer
    height = np.linalg.norm(axis_point - base_point) + 10.0

    print(
        f"Cylinder: base={base_point}, axis={axis_point}, radius={radius}, height={height}"
    )

    # Load structure and get residues
    residues = ribosome_entities(
        rcsb_id=rcsb_id,
        cifpath=cif_path,
        level="R",  # Residue level
        skip_nascent_chain=skip_nascent_chain,
    )

    print(f"Loaded {len(residues)} residues from structure")

    # Filter residues that fall within cylinder
    cylinder_residues = filter_residues_parallel(
        residues=residues,
        base_point=base_point,
        axis_point=axis_point,
        radius=radius,
        height=120,
    )

    print(f"Found {len(cylinder_residues)} residues within cylinder")

    # Group by chain and extract residue numbers
    cylinder_data = {}
    for residue in cylinder_residues:
        chain_id = residue.get_parent().id  # auth_asym_id
        residue_num = residue.id[1]  # auth_seq_id

        if chain_id not in cylinder_data:
            cylinder_data[chain_id] = []
        cylinder_data[chain_id].append(residue_num)

    # Sort residue numbers within each chain
    for chain_id in cylinder_data:
        cylinder_data[chain_id].sort()

    print(
        f"Cylinder residues by chain: {dict((k, len(v)) for k, v in cylinder_data.items())}"
    )

    # Write to JSON file
    with open(output_path, "w") as f:
        json.dump(cylinder_data, f, indent=2)

    print(f"Wrote cylinder residues to {output_path}")
    return cylinder_data



def create_npet_mesh(
    rcsb_id: str, log_dir: Optional[Path] = None, force: bool = False
) -> NPETProcessingTracker:
    """
    Creates NPET mesh for a given RCSB ID with comprehensive tracking and logging.

    This function orchestrates the execution of all pipeline stages, handling setup,
    error management, and tracking of the overall process.

    Args:
        rcsb_id: The RCSB PDB identifier
        log_dir: Directory to store logs (default: logs)
        force: Whether to force regeneration of existing assets

    Returns:
        NPETProcessingTracker: Processing tracker with detailed status information
    """
    print(f"Creating NPET mesh for {rcsb_id}")
    meshpath = AssetType.NPET_MESH.get_path(rcsb_id)
    if meshpath.exists():
        print(f"NPET mesh already exists for {rcsb_id}")

    rcsb_id = rcsb_id.upper()
    tracker = NPETProcessingTracker(rcsb_id, log_dir)
    
    try:
        # Initialize pipeline context to share data between stages
        context: Dict[str, Any] = {}
        
        # Define stage parameters
        alpha_shape_params = {
            "d3d_alpha": 200,
            "d3d_tol": 10,
            "d3d_offset": 3,

            "kdtree_radius": 40,
            "max_nn": 60,
            "tangent_planes_k": 20,

            "PR_depth": 6,
            "PR_ptweight": 4,
        }
        
        point_cloud_params = {
            "radius": 35,
            "height": 120,
            "voxel_size": 1,
            "atom_size": 2,
        }
        
        clustering_params = {
            "epsilon": 5.5,
            "min_samples": 600,
        }
        
        refinement_params = {
            "epsilon": 3.5,
            "min_samples": 175,
        }
        
        surface_params = {
            "alpha": 2,
            "tolerance": 1,
            "offset": 2,
        }
        
        normal_params = {
            "kdtree_radius": 10,
            "kdtree_max_nn": 15,
            "correction_tangent_planes_n": 10,
        }
        
        mesh_params = {
            "depth": 6,
            "ptweight": 3,
        }
        
        # Run the Setup stage
        setup_stage = SetupStage(rcsb_id, tracker, None, force)
        setup_result = setup_stage.run(context)
        context.update(setup_result)
        
        # Create artifacts directory
        artifacts_dir = Path(context["ro"].assets.paths.dir) / "artifacts"
        artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        # Run the PTC identification stage
        ptc_stage = PTCIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        ptc_result = ptc_stage.run(context)
        context.update(ptc_result)
        
        # Run the constriction identification stage
        constriction_stage = ConstrictionIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        constriction_result = constriction_stage.run(context)
        context.update(constriction_result)
        
        # Run the alpha shape stage
        alpha_stage = AlphaShapeStage(rcsb_id, tracker, artifacts_dir, True, alpha_shape_params)
        alpha_result = alpha_stage.run(context)
        context.update(alpha_result)
        
        # Run the landmark identification stage (verifies PTC and constriction)
        landmark_stage = LandmarkIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        landmark_result = landmark_stage.run(context)
        context.update(landmark_result)
        
        # Run the entity filtering stage
        entity_filtering_stage = EntityFilteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"]
        )
        entity_result = entity_filtering_stage.run(context)
        context.update(entity_result)
        
        # Run the point cloud processing stage
        point_cloud_stage = PointCloudProcessingStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"],
            voxel_size=point_cloud_params["voxel_size"],
            atom_size=point_cloud_params["atom_size"]
        )
        point_cloud_result = point_cloud_stage.run(context)
        context.update(point_cloud_result)
        
        # Run the clustering stage
        clustering_stage = ClusteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=clustering_params["epsilon"],
            min_samples=clustering_params["min_samples"]
        )
        clustering_result = clustering_stage.run(context)
        context.update(clustering_result)
        
        # Run the refinement stage
        refinement_stage = RefinementStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=refinement_params["epsilon"],
            min_samples=refinement_params["min_samples"]
        )
        refinement_result = refinement_stage.run(context)
        context.update(refinement_result)
        
        # Run the surface extraction stage
        surface_stage = SurfaceExtractionStage(
            rcsb_id, tracker, artifacts_dir, force,
            alpha=surface_params["alpha"],
            tolerance=surface_params["tolerance"],
            offset=surface_params["offset"]
        )
        surface_result = surface_stage.run(context)
        context.update(surface_result)
        
        # Run the normal estimation stage
        normal_stage = NormalEstimationStage(
            rcsb_id, tracker, artifacts_dir, force,
            kdtree_radius=normal_params["kdtree_radius"],
            kdtree_max_nn=normal_params["kdtree_max_nn"],
            correction_tangent_planes_n=normal_params["correction_tangent_planes_n"]
        )
        normal_result = normal_stage.run(context)
        context.update(normal_result)
        
        # Run the mesh reconstruction stage
        mesh_stage = MeshReconstructionStage(
            rcsb_id, tracker, artifacts_dir, force,
            depth=mesh_params["depth"],
            ptweight=mesh_params["ptweight"]
        )
        mesh_result = mesh_stage.run(context)
        context.update(mesh_result)
        
        # Run the validation stage
        validation_stage = ValidationStage(rcsb_id, tracker, artifacts_dir, force)
        validation_result = validation_stage.run(context)
        context.update(validation_result)
        
        # Complete pipeline and set final status
        tracker.complete_processing(True, watertight=context.get("watertight", False))
        
    except Exception as e:
        # Handle any uncaught exceptions
        if tracker.current_stage:
            tracker.end_stage(tracker.current_stage, False, error=e)
        tracker.complete_processing(False)
    
    return tracker