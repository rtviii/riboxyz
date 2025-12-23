import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import (
    NPETProcessingTracker,
    ProcessingStage,
)
from ribctl.lib.npet.various_visualization import (
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
)


def generate_clusters_pdb(
    clusters_container: dict, output_path: str, rcsb_id: str
) -> str:
    """
    Generate a PDB file with each DBSCAN cluster as a separate chain.
    """

    def get_chain_id(cluster_idx: int) -> str:
        """Generate unique chain IDs: A, B, C... Z, AA, AB, AC..."""
        if cluster_idx < 26:
            return chr(65 + cluster_idx)  # A-Z
        else:
            first = chr(65 + (cluster_idx - 26) // 26)
            second = chr(65 + (cluster_idx - 26) % 26)
            return first + second

    with open(output_path, "w") as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    DBSCAN CLUSTERS                         {rcsb_id}\n")
        pdb_file.write(f"TITLE     DBSCAN CLUSTERS FOR {rcsb_id}\n")
        pdb_file.write("REMARK    Generated from DBSCAN clustering\n")

        atom_num = 1
        valid_cluster_idx = 0
        connect_records = []  # Store CONECT records

        # Process each cluster (skip noise cluster -1)
        for cluster_id, points in clusters_container.items():
            if cluster_id == -1:  # Skip noise
                continue

            chain_id = get_chain_id(valid_cluster_idx)
            points_array = np.array(points)

            print(
                f"Writing cluster {cluster_id} as chain {chain_id} with {len(points_array)} points"
            )

            cluster_start_atom = atom_num

            # Write each point as a CA atom in the SAME residue (residue 1)
            for point_idx, (x, y, z) in enumerate(points_array):
                pdb_file.write(
                    f"ATOM  {atom_num:5d}  CA  ALA {chain_id}   1    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
                )

                # Create bonds between consecutive atoms in the cluster
                if point_idx > 0:
                    connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")

                atom_num += 1

            # Write chain terminator
            pdb_file.write(f"TER   {atom_num:5d}      ALA {chain_id}   1\n")
            atom_num += 1
            valid_cluster_idx += 1

        # Write all CONECT records at the end
        for connect in connect_records:
            pdb_file.write(connect)

        # Write end record
        pdb_file.write("END\n")

    print(f"Generated clusters PDB: {output_path}")
    print(f"Total clusters written: {valid_cluster_idx}")
    return output_path


class ClusteringStage(NPETPipelineStage):
    """
    Stage for applying DBSCAN clustering to identify the tunnel points.
    """

    def __init__(
        self,
        rcsb_id: str,
        tracker: NPETProcessingTracker,
        artifacts_dir: Path,
        force: bool = False,
        epsilon: float = 5.5,
        min_samples: int = 600,
    ):
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
        largest_cluster_path = (
            self.artifacts_dir / f"{self.rcsb_id}_largest_cluster.npy"
        )
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
                output_path=str(cluster_viz_path),
            )
            self.tracker.add_artifact(self.stage, cluster_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save cluster visualization: {str(viz_error)}")

        # Save individual clusters
        try:
            # Save each cluster as a separate file
            for cluster_id, cluster_points in clusters_container.items():
                if cluster_id != -1:  # Skip noise
                    cluster_path = (
                        self.artifacts_dir / f"{self.rcsb_id}_cluster_{cluster_id}.npy"
                    )
                    np.save(cluster_path, np.array(cluster_points))
                    if cluster_id == largest_cluster_id:
                        self.tracker.add_artifact(self.stage, cluster_path)
        except Exception as cluster_error:
            print(f"Warning: Could not save individual clusters: {str(cluster_error)}")

        return {
            "largest_cluster": largest_cluster,
            "clusters_container": clusters_container,
            "largest_cluster_id": largest_cluster_id,
        }
