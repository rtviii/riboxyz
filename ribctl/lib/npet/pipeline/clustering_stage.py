import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs


class ClusteringStage(NPETPipelineStage):
    """
    Stage for applying DBSCAN clustering to identify the tunnel points.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 epsilon: float = 5.5,
                 min_samples: int = 600):
        """
        Initialize the clustering stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            epsilon: DBSCAN epsilon parameter (neighborhood radius)
            min_samples: DBSCAN min_samples parameter (min points in neighborhood)
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.epsilon = epsilon
        self.min_samples = min_samples
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.CLUSTERING
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "epsilon": self.epsilon,
            "min_samples": self.min_samples,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Apply DBSCAN clustering to identify the tunnel points.
        
        Steps:
        1. Apply DBSCAN to the interior points
        2. Extract the largest cluster (representing the tunnel)
        3. Save the cluster data and visualization
        """
        empty_in_world_coords = context["empty_in_world_coords"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        
        # Apply DBSCAN clustering
        db, clusters_container = DBSCAN_capture(
            empty_in_world_coords, self.epsilon, self.min_samples
        )
        
        # Extract the largest cluster
        largest_cluster, largest_cluster_id = DBSCAN_pick_largest_cluster(
            clusters_container
        )
        
        # Save largest cluster
        largest_cluster_path = self.artifacts_dir / f"{self.rcsb_id}_largest_cluster.npy"
        np.save(largest_cluster_path, largest_cluster)
        self.tracker.add_artifact(self.stage, largest_cluster_path)
        
        # Save visualization if possible
        try:
            cluster_viz_path = self.artifacts_dir / f"{self.rcsb_id}_clusters.png"
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                clusters_container,
                self.epsilon,
                self.min_samples,
                ptc_pt,
                constriction_pt,
                largest_cluster,
                context.get("radius", 35),  # Default radius if not in context
                context.get("height", 120),  # Default height if not in context
                output_path=str(cluster_viz_path)
            )
            self.tracker.add_artifact(self.stage, cluster_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save cluster visualization: {str(viz_error)}")
        
        # Save individual clusters
        try:
            # Save each cluster as a separate file
            for cluster_id, cluster_points in clusters_container.items():
                if cluster_id != -1:  # Skip noise
                    cluster_path = self.artifacts_dir / f"{self.rcsb_id}_cluster_{cluster_id}.npy"
                    np.save(cluster_path, np.array(cluster_points))
                    if cluster_id == largest_cluster_id:
                        self.tracker.add_artifact(self.stage, cluster_path)
        except Exception as cluster_error:
            print(f"Warning: Could not save individual clusters: {str(cluster_error)}")
        
        return {
            "largest_cluster": largest_cluster,
            "clusters_container": clusters_container,
            "largest_cluster_id": largest_cluster_id
        }