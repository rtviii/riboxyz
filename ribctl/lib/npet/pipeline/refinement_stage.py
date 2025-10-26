import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.visualization.various_visualization import visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs


class RefinementStage(NPETPipelineStage):
    """
    Stage for refining the cluster to improve tunnel representation.
    
    This stage applies a second round of DBSCAN with stricter parameters
    to the largest cluster from the previous stage, in order to further 
    refine the tunnel representation.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 epsilon: float = 3.5,
                 min_samples: int = 175):
        """
        Initialize the refinement stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            epsilon: DBSCAN epsilon parameter for refinement
            min_samples: DBSCAN min_samples parameter for refinement
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.epsilon = epsilon
        self.min_samples = min_samples
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.REFINEMENT
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "epsilon": self.epsilon,
            "min_samples": self.min_samples,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply a second round of DBSCAN to refine the tunnel representation.
        
        Steps:
        1. Apply DBSCAN with stricter parameters to the largest cluster
        2. Extract the largest refined cluster
        3. Save the refined cluster data and visualization
        """
        largest_cluster = context["largest_cluster"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        
        # Apply second DBSCAN for refinement
        db_2, refined_clusters_container = DBSCAN_capture(
            largest_cluster, self.epsilon, self.min_samples
        )
        
        # Extract the largest refined cluster
        refined_cluster, refined_cluster_id = DBSCAN_pick_largest_cluster(
            refined_clusters_container
        )
        
        # Save refined cluster
        refined_cluster_path = self.artifacts_dir / f"{self.rcsb_id}_refined_cluster.npy"
        np.save(refined_cluster_path, refined_cluster)
        self.tracker.add_artifact(self.stage, refined_cluster_path)
        
        # Save visualization if possible
        try:
            refined_viz_path = self.artifacts_dir / f"{self.rcsb_id}_refined_clusters.png"
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                refined_clusters_container,
                self.epsilon,
                self.min_samples,
                ptc_pt,
                constriction_pt,
                refined_cluster,
                context.get("radius", 35),  # Default radius if not in context
                context.get("height", 120),  # Default height if not in context
                output_path=str(refined_viz_path)
            )
            self.tracker.add_artifact(self.stage, refined_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save refined cluster visualization: {str(viz_error)}")
        
        # Save individual refined clusters
        try:
            for cluster_id, cluster_points in refined_clusters_container.items():
                if cluster_id != -1:  # Skip noise
                    cluster_path = self.artifacts_dir / f"{self.rcsb_id}_refined_cluster_{cluster_id}.npy"
                    np.save(cluster_path, np.array(cluster_points))
                    if cluster_id == refined_cluster_id:
                        self.tracker.add_artifact(self.stage, cluster_path)
        except Exception as cluster_error:
            print(f"Warning: Could not save individual refined clusters: {str(cluster_error)}")
        
        return {
            "refined_cluster": refined_cluster,
            "refined_clusters_container": refined_clusters_container,
            "refined_cluster_id": refined_cluster_id
        }