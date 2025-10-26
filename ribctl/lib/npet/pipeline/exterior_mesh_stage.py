import os
import numpy as np
import pyvista as pv
import open3d as o3d
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.alphalib import cif_to_point_cloud, fast_normal_estimation, quick_surface_points, validate_mesh_pyvista
from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.visualization.various_visualization import visualize_pointcloud, visualize_mesh


class AlphaShapeStage(NPETPipelineStage):
    """
    Stage for generating the alpha shape representation of the ribosome structure.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 alpha_params: Dict[str, Any] = None):
        """
        Initialize the alpha shape stage with specific parameters.
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        if alpha_params is None:
            raise ValueError("Alpha shape parameters must be provided")
        self.params = alpha_params     
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.ALPHA_SHAPE
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return self.params
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Generate the alpha shape representation of the ribosome structure.
        """
        # Get paths from context
        cifpath = context["cifpath"]
        ashapepath = context["ashapepath"]
        
        # Check if alpha shape already exists
        regenerate_alpha = self.force or not os.path.exists(ashapepath)
        
        if not regenerate_alpha:
            print(f"Using existing alpha shape mesh: {ashapepath}")
            self.tracker.add_artifact(self.stage, Path(ashapepath))
        else:
            print(f"Generating alpha shape mesh for {self.rcsb_id}")
            
            # Generate point cloud
            ptcloudpath = self.artifacts_dir / f"{self.rcsb_id}_structure_ptcloud.npy"
            
            if not os.path.exists(ptcloudpath) or self.force:
                print("Extracting point cloud from CIF file")
                first_assembly_chains = context["ro"].first_assembly_auth_asym_ids()
                ptcloud = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
                np.save(ptcloudpath, ptcloud)
            else:
                print(f"Loading existing point cloud from {ptcloudpath}")
                ptcloud = np.load(ptcloudpath)
            
            self.tracker.add_artifact(self.stage, ptcloudpath)
            
            # Surface extraction
            print("Beginning Delaunay 3D reconstruction")
            surface_pts = quick_surface_points(
                ptcloud, 
                self.params["d3d_alpha"], 
                self.params["d3d_tol"], 
                self.params["d3d_offset"]
            )
            
            # Save surface points
            surface_pts_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_surface_points.npy"
            np.save(surface_pts_path, surface_pts)
            self.tracker.add_artifact(self.stage, surface_pts_path)
            
            # Visualize if possible
            try:
                surface_viz_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_surface.png"
                visualize_pointcloud(surface_pts, self.rcsb_id, output_path=str(surface_viz_path))
                self.tracker.add_artifact(self.stage, surface_viz_path)
            except Exception as viz_error:
                print(f"Warning: Could not save visualization: {str(viz_error)}")
            
            # Normal estimation
            normal_estimated_pcd = fast_normal_estimation(
                surface_pts, 
                self.params["kdtree_radius"], 
                self.params["max_nn"], 
                self.params["tangent_planes_k"]
            )
            
            # Save normal-estimated point cloud
            alpha_normals_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_normals.ply"
            o3d.io.write_point_cloud(str(alpha_normals_path), normal_estimated_pcd)
            self.tracker.add_artifact(self.stage, alpha_normals_path)
            
            # Apply Poisson reconstruction
            apply_poisson_reconstruction(
                str(alpha_normals_path),
                ashapepath,
                recon_depth=self.params["PR_depth"],
                recon_pt_weight=self.params["PR_ptweight"],
            )
            
            # Extract largest component for better quality
            mesh = pv.read(ashapepath)
            labeled = mesh.connectivity(largest=True)
            labeled.save(ashapepath)
            
            # Check watertightness
            watertight = validate_mesh_pyvista(labeled)
            
            if not watertight:
                print("Warning: Alpha shape mesh is not watertight, but continuing anyway")
            
            # Save mesh visualization
            try:
                mesh_viz_path = self.artifacts_dir / f"{self.rcsb_id}_alpha_mesh.png"
                visualize_mesh(ashapepath, self.rcsb_id, output_path=str(mesh_viz_path))
                self.tracker.add_artifact(self.stage, mesh_viz_path)
            except Exception as viz_error:
                print(f"Warning: Could not save visualization: {str(viz_error)}")
            
            # Also add ASCII version as artifact if it exists
            ascii_path = Path(str(ashapepath).split(".")[0] + "_ascii.ply")
            if ascii_path.exists():
                self.tracker.add_artifact(self.stage, ascii_path)
        
        # Add the final alpha shape as artifact
        self.tracker.add_artifact(self.stage, Path(ashapepath))
        
        # Verify alpha shape exists for the remaining pipeline
        if not os.path.exists(ashapepath):
            raise FileNotFoundError(f"Alpha shape file {ashapepath} not found after generation")
            
        return {
            "ashapepath": ashapepath
        }