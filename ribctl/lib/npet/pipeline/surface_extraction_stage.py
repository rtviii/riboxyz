import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import ptcloud_convex_hull_points
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.visualization.various_visualization import visualize_pointcloud


class SurfaceExtractionStage(NPETPipelineStage):
    """
    Stage for extracting the surface points from the refined cluster.
    
    This stage extracts points on the surface of the refined cluster
    using a convex hull approach, which will be used for mesh generation.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 alpha: float = 2,
                 tolerance: float = 1,
                 offset: float = 2):
        """
        Initialize the surface extraction stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            alpha: Alpha value for surface extraction
            tolerance: Tolerance for surface extraction
            offset: Offset for surface extraction
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.alpha = alpha
        self.tolerance = tolerance
        self.offset = offset
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.SURFACE_EXTRACTION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "alpha": self.alpha,
            "tolerance": self.tolerance,
            "offset": self.offset,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Extract surface points from the refined cluster.
        
        Steps:
        1. Apply convex hull algorithm to extract surface points
        2. Save the surface points and visualization
        """
        refined_cluster = context["refined_cluster"]
        
        # Extract surface points
        surface_pts = ptcloud_convex_hull_points(
            refined_cluster, self.alpha, self.tolerance, self.offset
        )
        
        # Save surface points
        surface_pts_path = self.artifacts_dir / f"{self.rcsb_id}_surface_points.npy"
        np.save(surface_pts_path, surface_pts)
        self.tracker.add_artifact(self.stage, surface_pts_path)
        
        # Save visualization if possible
        try:
            surface_viz_path = self.artifacts_dir / f"{self.rcsb_id}_surface_points.png"
            visualize_pointcloud(surface_pts, self.rcsb_id, output_path=str(surface_viz_path))
            self.tracker.add_artifact(self.stage, surface_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save surface point visualization: {str(viz_error)}")
        
        return {
            "surface_pts": surface_pts
        }