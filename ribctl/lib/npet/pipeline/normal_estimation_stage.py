import open3d as o3d
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import estimate_normals
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage


class NormalEstimationStage(NPETPipelineStage):
    """
    Stage for estimating normals for the surface points.
    
    This stage estimates surface normals for the point cloud,
    which are required for the Poisson surface reconstruction.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 kdtree_radius: float = 10,
                 kdtree_max_nn: int = 15,
                 correction_tangent_planes_n: int = 10):
        """
        Initialize the normal estimation stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            kdtree_radius: Radius for kdtree search during normal estimation
            kdtree_max_nn: Maximum number of neighbors for normal estimation
            correction_tangent_planes_n: Number of neighbors for tangent plane fitting
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.kdtree_radius = kdtree_radius
        self.kdtree_max_nn = kdtree_max_nn
        self.correction_tangent_planes_n = correction_tangent_planes_n
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.NORMAL_ESTIMATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "kdtree_radius": self.kdtree_radius,
            "kdtree_max_nn": self.kdtree_max_nn,
            "correction_tangent_planes_n": self.correction_tangent_planes_n,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Estimate normals for the surface points.
        
        Steps:
        1. Apply normal estimation to the surface points
        2. Save the point cloud with normals
        """
        surface_pts = context["surface_pts"]
        
        # Estimate normals
        normal_estimated_pcd = estimate_normals(
            surface_pts,
            kdtree_radius=self.kdtree_radius,
            kdtree_max_nn=self.kdtree_max_nn,
            correction_tangent_planes_n=self.correction_tangent_planes_n,
        )
        
        # Save the point cloud with normals
        normals_pcd_path = self.artifacts_dir / f"{self.rcsb_id}_normal_estimated_pcd.ply"
        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        self.tracker.add_artifact(self.stage, normals_pcd_path)
        
        return {
            "normal_estimated_pcd": normal_estimated_pcd,
            "normals_pcd_path": normals_pcd_path
        }