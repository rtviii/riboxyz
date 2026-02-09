# ribctl/lib/npet2/core/config.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Literal, Tuple


@dataclass(frozen=True)
class GridLevelConfig:
    name: str
    voxel_size_A: float
    atom_radius_mode: Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0
    occupancy_backend: Literal["legacy_kdtree", "edt"] = "legacy_kdtree"


@dataclass(frozen=True)
class RunConfig:
    # === Chain selection ===
    occupancy_chain_mode: Literal["walls_only", "assembly_all"] = "walls_only"
    occupancy_exclude_trna: bool = True
    occupancy_exclude_auth_asym_ids: Tuple[str, ...] = ()

    # === Region definition ===
    cylinder_radius_A: float = 35.0
    cylinder_height_A: float = 120.0

    # === Stage20: Exterior shell (whole ribosome surface) ===
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

    # === Stage50: DBSCAN clustering on level_0 (1.0A grid) ===
    dbscan_level0_coarse_eps_A: float = 5.5
    dbscan_level0_coarse_min_samples: int = 600
    dbscan_level0_refine_eps_A: float = 3.5
    dbscan_level0_refine_min_samples: int = 175
    mesh_level0_enable: bool = True

    # === Stage55: Grid refinement (0.5A ROI pass) ===
    refine_voxel_size_A: float = 0.5
    refine_roi_pad_A: float = 10.0
    refine_atom_radius_A: float = 2.0
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

    mesh_level1_enable: bool = True

    # === Meshing (MC + smoothing, shared by Stage50/55/70) ===
    mesh_gaussian_sigma_voxels: float = 1.5
    mesh_smooth_method        : str   = "taubin"
    mesh_taubin_pass_band     : float = 0.1
    mesh_level0_smooth_iters  : int   = 20
    mesh_level1_smooth_iters  : int   = 40
    mesh_fill_holes_A         : float = 100.0
    mesh_atom_clearance_A     : float = 1.5