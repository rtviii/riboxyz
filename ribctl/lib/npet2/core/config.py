# ribctl/lib/npet2/core/config.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional


@dataclass(frozen=True)
class DBSCANPolicy:
    eps_A: float                     # physical eps in Å
    density_factor: float = 0.10     # α in your back-of-envelope
    min_samples_override: Optional[int] = None


from dataclasses import dataclass, field
from typing import List, Literal


@dataclass(frozen=True)
class GridLevelConfig:
    name: str
    voxel_size_A: float

    atom_radius_mode: Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0

    occupancy_backend: Literal["legacy_kdtree", "grid_stamp", "edt", "gpu"] = "legacy_kdtree"
    roi_backend: Literal["full_cylinder", "bbox_from_prev", "tube_from_prev"] = "full_cylinder"


@dataclass(frozen=True)
class RunConfig:
    # region definition
    cylinder_radius_A: float = 35.0
    cylinder_height_A: float = 120.0

    # exterior shell (legacy alpha stage parameters)
    alpha_d3d_alpha: float = 200
    alpha_d3d_tol: float = 10
    alpha_d3d_offset: float = 3
    alpha_kdtree_radius: float = 40
    alpha_max_nn: int = 60
    alpha_tangent_planes_k: int = 20
    alpha_poisson_depth: int = 6
    alpha_poisson_ptweight: int = 4
    alpha_fill_holes: float = 2000

    # clustering (legacy values)
    dbscan_eps_A: float = 5.5
    dbscan_min_samples: int = 600
    refine_eps_A: float = 3.5
    refine_min_samples: int = 175

    # surface extraction
    surface_alpha: float = 2
    surface_tolerance: float = 1
    surface_offset: float = 2

    # normals
    normals_radius: float = 10
    normals_max_nn: int = 15
    normals_tangent_k: int = 10

    # mesh reconstruction
    mesh_poisson_depth: int = 6
    mesh_poisson_ptweight: int = 3

    # --- grid refinement (0.5Å ROI pass)
    refine_voxel_size_A      : float = 0.5
    refine_roi_pad_A         : float = 10.0
    refine_atom_radius_A     : float = 2.0
    refine_cc_connectivity   : int   = 26
    refine_topk_preview      : int   = 5
    refine_max_preview_points: int   = 50_000


    # refinement plan (kept for later; default only one level)
    grid_levels: List[GridLevelConfig] = field(default_factory=lambda: [
        GridLevelConfig(name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"),
    ])

