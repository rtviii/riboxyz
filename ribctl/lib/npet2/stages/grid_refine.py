# ribctl/lib/npet2/stages/grid_refine.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pyvista as pv
from scipy import ndimage
import time


from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    transform_points_from_C0,
)

from ribctl.lib.npet2.backends.grid_occupancy import (
    GridSpec,
    occupancy_via_edt,
    connected_components_3d,
)


def _make_bbox_grid(lo: np.ndarray, hi: np.ndarray, voxel: float) -> GridSpec:
    lo = np.asarray(lo, dtype=np.float32)
    hi = np.asarray(hi, dtype=np.float32)
    voxel = float(voxel)

    span = hi - lo
    # +1 so both ends are representable
    shape = tuple((np.ceil(span / voxel).astype(np.int32) + 1).tolist())
    return GridSpec(origin=lo, voxel_size=voxel, shape=shape)


def _voxel_centers_from_indices(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin[None, :] + ijk * float(grid.voxel_size)


def _points_to_ijk(grid: GridSpec, pts_c0: np.ndarray) -> np.ndarray:
    """Nearest-voxel mapping for points in C0 -> ijk indices."""
    v = float(grid.voxel_size)
    ijk = np.floor((pts_c0 - grid.origin[None, :]) / v + 0.5).astype(np.int32)
    return ijk


def _valid_ijk(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    nx, ny, nz = grid.shape
    m = (
        (ijk[:, 0] >= 0)
        & (ijk[:, 0] < nx)
        & (ijk[:, 1] >= 0)
        & (ijk[:, 1] < ny)
        & (ijk[:, 2] >= 0)
        & (ijk[:, 2] < nz)
    )
    return m


def _topk_component_stats(labeled: np.ndarray, k: int = 6) -> list[dict]:
    sizes = np.bincount(labeled.ravel())
    if sizes.size == 0:
        return []
    sizes[0] = 0
    if sizes.sum() == 0:
        return []

    top_labels = np.argsort(sizes)[::-1][:k]
    out = []
    for lab in top_labels:
        if lab == 0 or sizes[lab] == 0:
            continue
        idx = np.argwhere(labeled == lab)
        bbox_min = idx.min(axis=0)
        bbox_max = idx.max(axis=0)
        out.append(
            {
                "label": int(lab),
                "size": int(sizes[lab]),
                "bbox_min": bbox_min.tolist(),
                "bbox_max": bbox_max.tolist(),
            }
        )
    return out
