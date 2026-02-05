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

    # === Stage20: Exterior shell (legacy alpha shape) ===
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

    surface_alpha: float = 200
    surface_tolerance: float = 10
    surface_offset: float = 3
    
    normals_radius: float = 10
    normals_max_nn: int = 15
    normals_tangent_k: int = 10
    # === Stage50: DBSCAN clustering on level_0 (1.0Å grid) ===
    # Coarse pass: merge large regions, bridge small gaps
    # Rule: eps ≈ 2-3 × voxel_size, min_samples ≈ 30-50
    dbscan_level0_coarse_eps_A: float = 5.5  # ~5.5× for legacy compatibility
    dbscan_level0_coarse_min_samples: int = 600  # high for robustness

    # Refine pass: tighten on largest cluster from coarse
    # Rule: eps ≈ 1.5-2 × voxel_size (relative to coarse cluster density)
    dbscan_level0_refine_eps_A: float = 3.5
    dbscan_level0_refine_min_samples: int = 175

    # Mesh generation (optional, for Stage50)
    mesh_level0_enable: bool = True
    mesh_level0_poisson_depth: int = 6
    mesh_level0_poisson_ptweight: int = 3

    # === Stage55: Grid refinement (0.5Å ROI pass) ===
    refine_voxel_size_A: float = 0.5
    refine_roi_pad_A: float = 10.0
    refine_atom_radius_A: float = 2.5  # slightly larger to seal cracks

    # Localization: restrict refined void to vicinity of Stage50 winner
    refine_keep_within_A: float = 6.0  # distance from coarse tunnel (0=off)

    # Morphology cleanup (before DBSCAN)
    refine_occ_close_iters: int = 0  # seal occupancy cracks
    refine_void_open_iters: int = 1  # break skinny bridges in void

    # Prevent ROI boundary contamination
    refine_forbid_roi_boundary: bool = True

    # DBSCAN on refined boundary points (0.5Å grid)
    # Coarse pass: eps ≈ 6× voxel_size (legacy compat from original Stage55)
    dbscan_level1_coarse_eps_A: float = 3.0
    dbscan_level1_coarse_min_samples: int = 30

    # Refine pass: tighten further
    dbscan_level1_refine_eps_A: float = 3.0
    dbscan_level1_refine_min_samples: int = 20

    # Safety caps for DBSCAN performance
    refine_dbscan_max_points: int = 0  # subsample if >0
    refine_dbscan_seed: int = 0

    # Diagnostics
    refine_dbscan_max_cluster_stats: int = 25

    # Mesh generation (Stage55)
    mesh_level1_enable: bool = True

    # === Stage60: Surface normals ===
    normals_radius: float = 10
    normals_max_nn: int = 15
    normals_tangent_k: int = 10

    # === Stage70: Final mesh validation ===
    # Voxel-based fallback mesh cleanup
    voxel_mesh_fill_holes_A: float = 50.0
    voxel_mesh_smooth_iters: int = 0

    mesh_poisson_depth: int = 6
    mesh_poisson_ptweight: int = 3
    
    # Voxel-based fallback mesh cleanup
    voxel_mesh_fill_holes_A: float = 50.0
    voxel_mesh_smooth_iters: int = 0

    mesh_level1_poisson_depth: int = 8       # ← deeper for detail
    mesh_level1_poisson_ptweight: float = 0.5  # ← LOWER for smoothness
    
    # Stage70 voxel mesh smoothing
    voxel_mesh_smooth_iters: int = 10  # ← increase from 0 to 10 for smooth simulation mesh