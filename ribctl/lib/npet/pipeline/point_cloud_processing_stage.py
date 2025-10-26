import numpy as np
import pyvista as pv
from pathlib import Path
from typing import Any, Dict

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    create_point_cloud_mask,
    transform_points_from_C0,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import (
    NPETProcessingTracker,
    ProcessingStage,
)
from ribctl.lib.npet.visualization.various_visualization import visualize_pointcloud


class PointCloudProcessingStage(NPETPipelineStage):
    """
    Stage for transforming points to C0 space, creating the tunnel mask,
    and generating points representing the empty space inside the ribosome.
    """

    def __init__(
        self,
        rcsb_id: str,
        tracker: NPETProcessingTracker,
        artifacts_dir: Path,
        force: bool = False,
        radius: float = 35,
        height: float = 120,
        voxel_size: float = 1,
        atom_size: float = 2,
    ):
        """
        Initialize the point cloud processing stage.

        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            radius: Radius of the tunnel cylinder
            height: Height of the tunnel cylinder
            voxel_size: Size of voxels for discretization
            atom_size: Size of atoms for creating the mask
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.radius = radius
        self.height = height
        self.voxel_size = voxel_size
        self.atom_size = atom_size

    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.POINT_CLOUD_PROCESSING

    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "radius": self.radius,
            "height": self.height,
            "voxel_size": self.voxel_size,
            "atom_size": self.atom_size,
        }

    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Process the filtered points to create a representation of the tunnel.

        Steps:
        1. Transform filtered points to canonical space (C0)
        2. Create a mask representing occupied space
        3. Find points in empty space within the tunnel region
        4. Transform back to world coordinates
        5. Select only points inside the alpha shape mesh
        """
        filtered_points = context["filtered_points"]
        ptc_pt           = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        ashapepath = context["ashapepath"]

        # Transform points to canonical space (C0)
        transformed_points = transform_points_to_C0(
            filtered_points, ptc_pt, constriction_pt
        )

        # Save transformed points
        transformed_points_path = (
            self.artifacts_dir / f"{self.rcsb_id}_transformed_points.npy"
        )
        np.save(transformed_points_path, transformed_points)
        self.tracker.add_artifact(self.stage, transformed_points_path)

        # Create mask representing occupied space
        mask, (x, y, z) = create_point_cloud_mask(
            transformed_points,
            radius=self.radius,
            height=self.height,
            voxel_size=self.voxel_size,
            radius_around_point=self.atom_size,
        )

        # Find coordinates of empty space
        points = np.where(~mask)
        empty_coordinates = np.column_stack((x[points[0]], y[points[1]], z[points[2]]))

        # Transform back to world coordinates
        back_projected = transform_points_from_C0(
            empty_coordinates, ptc_pt, constriction_pt
        )

        # Save back projected points
        back_projected_path = self.artifacts_dir / f"{self.rcsb_id}_back_projected.npy"
        np.save(back_projected_path, back_projected)
        self.tracker.add_artifact(self.stage, back_projected_path)

        ashape_watertight_mesh = pv.read(ashapepath)
        select = pv.PolyData(back_projected).select_enclosed_points(
            ashape_watertight_mesh
        )
        mask = select["SelectedPoints"]
        interior = back_projected[mask == 1]
        empty_in_world_coords = np.array(interior)

        # Save interior points
        interior_points_path = (
            self.artifacts_dir / f"{self.rcsb_id}_interior_points.npy"
        )
        np.save(interior_points_path, empty_in_world_coords)
        self.tracker.add_artifact(self.stage, interior_points_path)

        # Save visualization if possible
        try:
            pc_viz_path = self.artifacts_dir / f"{self.rcsb_id}_point_cloud.png"
            visualize_pointcloud(
                empty_in_world_coords, self.rcsb_id, output_path=str(pc_viz_path)
            )
            self.tracker.add_artifact(self.stage, pc_viz_path)
        except Exception as viz_error:
            print(
                f"Warning: Could not save point cloud visualization: {str(viz_error)}"
            )

        return {"empty_in_world_coords": empty_in_world_coords}
