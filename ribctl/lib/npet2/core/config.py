# ribctl/lib/npet2/core/config.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional


@dataclass(frozen=True)
class GridLevelConfig:
    name: str
    voxel_size_A: float

    atom_radius_mode: Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0

    occupancy_backend: Literal["legacy_kdtree", "grid_stamp", "edt", "gpu"] = (
        "legacy_kdtree"
    )
    roi_backend: Literal["full_cylinder", "bbox_from_prev", "tube_from_prev"] = (
        "full_cylinder"
    )


@dataclass(frozen=True)
class RunConfig:
    # === Region definition ===
    cylinder_radius_A: float = 35.0
    cylinder_height_A: float = 120.0

    # === Stage20: Exterior shell (whole ribosome surface) ===
    # Large alpha because the ribosome is roughly convex
    alpha_d3d_alpha: float = 200
    alpha_d3d_tol: float = 10
    alpha_d3d_offset: float = 3
    alpha_kdtree_radius: float = 40
    alpha_max_nn: int = 60
    alpha_tangent_planes_k: int = 20
    alpha_poisson_depth: int = 6
    alpha_poisson_ptweight: int = 4
    alpha_fill_holes: float = 2000

    # === Stage40: Grid levels ===
    grid_levels: List[GridLevelConfig] = field(
        default_factory=lambda: [
            GridLevelConfig(
                name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"
            ),
        ]
    )

    # === Tunnel surface extraction (Delaunay on DBSCAN cluster) ===
    # Small alpha because the tunnel is deeply concave
    tunnel_surface_alpha: float = 2.0
    tunnel_surface_tolerance: float = 1.0
    tunnel_surface_offset: float = 2.0

    # === Surface normals (shared by Stage50 mesh, Stage60, Stage55) ===
    normals_radius: float = 10
    normals_max_nn: int = 15
    normals_tangent_k: int = 10

    # === Stage50: DBSCAN clustering on level_0 (1.0A grid) ===
    dbscan_level0_coarse_eps_A: float = 5.5
    dbscan_level0_coarse_min_samples: int = 600
    dbscan_level0_refine_eps_A: float = 3.5
    dbscan_level0_refine_min_samples: int = 175

    # Stage50 mesh generation
    mesh_level0_enable: bool = True
    mesh_level0_poisson_depth: int = 6
    mesh_level0_poisson_ptweight: float = 3.0

    # === Stage55: Grid refinement (0.5A ROI pass) ===
    refine_voxel_size_A: float = 0.5
    refine_roi_pad_A: float = 10.0
    refine_atom_radius_A: float = 2.5

    refine_keep_within_A: float = 6.0
    refine_occ_close_iters: int = 0
    refine_void_open_iters: int = 1
    refine_forbid_roi_boundary: bool = True

    # DBSCAN on refined grid
    dbscan_level1_coarse_eps_A: float = 3.0
    dbscan_level1_coarse_min_samples: int = 30
    dbscan_level1_refine_eps_A: float = 3.0
    dbscan_level1_refine_min_samples: int = 20

    refine_dbscan_max_points: int = 0
    refine_dbscan_seed: int = 0
    refine_dbscan_max_cluster_stats: int = 25

    # Stage55 mesh generation
    mesh_level1_enable: bool = True
    mesh_level1_poisson_depth: int = 8
    mesh_level1_poisson_ptweight: float = 3.0

    # === Stage70: Final mesh validation ===
    mesh_poisson_depth: int = 6
    mesh_poisson_ptweight: float = 3.0

    voxel_mesh_fill_holes_A: float = 50.0
    voxel_mesh_smooth_iters: int = 10
