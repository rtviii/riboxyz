# npet2/core/config.py
"""
All npet2 configuration in one place.

- Settings: paths, env vars, external tool locations
- GridLevelConfig / RunConfig: numerical pipeline parameters

Everything is overridable via environment variables (for Docker)
or programmatically.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Literal, Tuple



def _env_path(var: str, default: str) -> Path:
    return Path(os.environ.get(var, default))


@dataclass(frozen=True)
class Settings:
    npet2_root:        Path = field(default_factory=lambda: _env_path("NPET2_ROOT", str(Path.home() / "npet2_data")))
    runs_root:         Path = field(default_factory=lambda: _env_path("NPET2_RUNS_ROOT", str(_env_path("NPET2_ROOT", str(Path.home() / "npet2_data")) / "runs")))
    cache_root:        Path = field(default_factory=lambda: _env_path("NPET2_CACHE_ROOT", str(_env_path("NPET2_ROOT", str(Path.home() / "npet2_data")) / "cache")))
    poisson_recon_bin: str  = field(default_factory=lambda: os.environ.get("NPET2_POISSON_RECON_BIN", "PoissonRecon"))
    riboxyz_api_base:  str  = field(default_factory=lambda: os.environ.get("NPET2_RIBOXYZ_API_URL", "http://localhost:8000"))


# Module-level singleton. Import this wherever you need paths.
SETTINGS = Settings()


# ---------------------------------------------------------------------------
# Grid level config
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class GridLevelConfig:
    name                 : str
    voxel_size_A         : float
    atom_radius_mode     : Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0
    occupancy_backend    : Literal["legacy_kdtree", "edt"] = "legacy_kdtree"


# ---------------------------------------------------------------------------
# Run config (all numerical pipeline parameters)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RunConfig:
    # === Chain selection ===
    occupancy_chain_mode: Literal["walls_only", "assembly_all"] = "walls_only"
    occupancy_exclude_trna: bool = True
    occupancy_exclude_auth_asym_ids: Tuple[str, ...] = ()

    # === Region definition ===
    cylinder_radius_A: float = 35
    cylinder_height_A: float = 120
    cylinder_ptc_extension_A: float = 20

    # === Stage20: Exterior shell (whole ribosome surface) ===
    alpha_d3d_alpha       : float = 200
    alpha_d3d_tol         : float = 10
    alpha_d3d_offset      : float = 3
    alpha_kdtree_radius   : float = 40
    alpha_max_nn          : int   = 60
    alpha_tangent_planes_k: int   = 20
    alpha_poisson_depth   : int   = 6
    alpha_poisson_ptweight: int   = 4
    alpha_fill_holes      : float = 2000

    # === Stage40: Grid levels ===
    grid_levels: List[GridLevelConfig] = field(
        default_factory=lambda: [
            GridLevelConfig(
                name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"
            ),
        ]
    )

    # === Stage50: DBSCAN clustering on level_0 ===
    dbscan_level0_coarse_eps_A      : float = 5.5
    dbscan_level0_coarse_min_samples: int   = 600
    dbscan_level0_refine_eps_A      : float = 3.5
    dbscan_level0_refine_min_samples: int   = 175
    mesh_level0_enable              : bool  = True

    # === Stage55: Grid refinement (0.5A ROI pass) ===
    refine_voxel_size_A       : float = 0.5
    refine_roi_pad_A          : float = 10.0
    refine_atom_radius_A      : float = 2.0
    refine_keep_within_A      : float = 6.0
    refine_occ_close_iters    : int   = 0
    refine_void_open_iters    : int   = 1
    refine_forbid_roi_boundary: bool  = True

    dbscan_level1_coarse_eps_A      : float = 3.0
    dbscan_level1_coarse_min_samples: int   = 30
    dbscan_level1_refine_eps_A      : float = 3.0
    dbscan_level1_refine_min_samples: int   = 20

    refine_dbscan_max_points        : int   = 0
    refine_dbscan_seed              : int   = 0

    mesh_level1_enable: bool = True

    # === Meshing (MC + smoothing) ===
    mesh_smooth_method        : str   = "taubin"
    mesh_level0_gaussian_sigma: float = 1.0
    mesh_level1_gaussian_sigma: float = 1.5
    mesh_taubin_pass_band     : float = 0.1

    mesh_level0_smooth_iters  : int   = 40
    mesh_level1_smooth_iters  : int   = 60

    mesh_fill_holes_A         : float = 100.0
    mesh_atom_clearance_A     : float = 1.5