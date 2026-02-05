from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class GridSpec:
    # origin and voxel size define world coords: x = origin + i*voxel
    origin: np.ndarray          # (3,)
    voxel_size: float
    shape: Tuple[int, int, int] # (nx, ny, nz)


def make_cylinder_grid(radius_A: float, height_A: float, voxel_A: float) -> GridSpec:
    """
    Canonical cylinder in C0:
      x in [-R, R], y in [-R, R], z in [0, H]
    """
    nx = int(np.floor((2 * radius_A) / voxel_A)) + 1
    ny = int(np.floor((2 * radius_A) / voxel_A)) + 1
    nz = int(np.floor(height_A / voxel_A)) + 1
    origin = np.array([-radius_A, -radius_A, 0.0], dtype=np.float32)
    return GridSpec(origin=origin, voxel_size=float(voxel_A), shape=(nx, ny, nz))


def grid_world_coords(grid: GridSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns axis coordinate arrays (x, y, z) for voxel centers along each axis.
    """
    ox, oy, oz = grid.origin
    nx, ny, nz = grid.shape
    v = grid.voxel_size
    x = ox + np.arange(nx, dtype=np.float32) * v
    y = oy + np.arange(ny, dtype=np.float32) * v
    z = oz + np.arange(nz, dtype=np.float32) * v
    return x, y, z


def cylinder_mask(grid: GridSpec, radius_A: float) -> np.ndarray:
    """
    Boolean mask of voxels inside cylinder radius (in C0).
    """
    x, y, _ = grid_world_coords(grid)
    X, Y = np.meshgrid(x, y, indexing="ij")
    inside = (X * X + Y * Y) <= (radius_A * radius_A)
    # broadcast across z
    return inside[:, :, None]


def points_to_occupied_seeds(points_c0: np.ndarray, grid: GridSpec) -> np.ndarray:
    """
    Convert points (C0 coords) to a sparse occupied grid of seeds at nearest voxels.
    """
    pts = np.asarray(points_c0, dtype=np.float32)
    v = grid.voxel_size
    origin = grid.origin

    ijk = np.floor((pts - origin[None, :]) / v + 0.5).astype(np.int32)
    nx, ny, nz = grid.shape

    valid = (
        (ijk[:, 0] >= 0) & (ijk[:, 0] < nx) &
        (ijk[:, 1] >= 0) & (ijk[:, 1] < ny) &
        (ijk[:, 2] >= 0) & (ijk[:, 2] < nz)
    )
    ijk = ijk[valid]
    occ = np.zeros(grid.shape, dtype=np.bool_)
    if ijk.shape[0] > 0:
        occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return occ


def occupancy_via_edt(points_c0: np.ndarray, grid: GridSpec, atom_radius_A: float) -> np.ndarray:
    """
    Occupancy grid: voxel is occupied if within atom_radius_A of any atom center.

    Steps:
      - seed occupied at nearest voxels
      - edt on ~occupied gives distance (in voxels) to nearest seed
      - threshold <= r_vox
    """
    seeds = points_to_occupied_seeds(points_c0, grid)

    # Distance (in voxels) from each voxel to nearest True in 'seeds'
    # distance_transform_edt computes distance to nearest zero;
    # so we compute on ~seeds, where zeros correspond to seeds.
    dist_vox = ndimage.distance_transform_edt(~seeds)

    r_vox = float(atom_radius_A) / float(grid.voxel_size)
    occupied = dist_vox <= r_vox
    return occupied


def empty_points_from_mask(grid: GridSpec, empty_mask: np.ndarray) -> np.ndarray:
    """
    Return coordinates (C0) of voxel centers where empty_mask is True.
    """
    empty_idx = np.where(empty_mask)
    x, y, z = grid_world_coords(grid)
    pts = np.column_stack((x[empty_idx[0]], y[empty_idx[1]], z[empty_idx[2]])).astype(np.float32)
    return pts


def save_grid_npy(grid: GridSpec, data: np.ndarray, path: Path, *, compress: bool = False) -> None:
    """
    Save a 3D grid along with its GridSpec metadata.
    File format: {path}_data.npy + {path}_spec.json
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if compress:
        np.savez_compressed(data_path.with_suffix('.npz'), data=data)
    else:
        np.save(data_path, data)
    
    spec_dict = {
        "origin": grid.origin.tolist(),
        "voxel_size": float(grid.voxel_size),
        "shape": list(grid.shape),
    }
    spec_path.write_text(__import__('json').dumps(spec_dict, indent=2))


def load_grid_npy(path: Path) -> tuple[GridSpec, np.ndarray]:
    """
    Load a grid saved by save_grid_npy.
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if data_path.with_suffix('.npz').exists():
        data = np.load(data_path.with_suffix('.npz'))['data']
    else:
        data = np.load(data_path)
    
    spec_dict = __import__('json').loads(spec_path.read_text())
    grid = GridSpec(
        origin=np.array(spec_dict["origin"], dtype=np.float32),
        voxel_size=float(spec_dict["voxel_size"]),
        shape=tuple(spec_dict["shape"]),
    )
    return grid, data


def voxel_to_world(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    """
    Convert voxel indices (i,j,k) to world coordinates.
    ijk: (N, 3) or (3,) array of voxel indices
    Returns: (N, 3) or (3,) world coordinates
    """
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin + ijk * grid.voxel_size


def world_to_voxel(grid: GridSpec, xyz: np.ndarray) -> np.ndarray:
    """
    Convert world coordinates to voxel indices.
    xyz: (N, 3) or (3,) world coordinates
    Returns: (N, 3) or (3,) voxel indices (floats; use floor/round as needed)
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    return (xyz - grid.origin) / grid.voxel_size


def get_occupied_voxel_centers(grid: GridSpec, occupancy: np.ndarray) -> np.ndarray:
    """
    Get world coordinates of occupied voxel centers.
    """
    occupied_idx = np.argwhere(occupancy)
    return voxel_to_world(grid, occupied_idx)
# ribctl/lib/npet2/backends/grid_occupancy.py
# ADD these functions at the end:

from scipy import ndimage

def connected_components_3d(
    binary_mask: np.ndarray, 
    connectivity: int = 26
) -> tuple[np.ndarray, int]:
    """
    Find connected components in a 3D binary mask.
    
    Args:
        binary_mask: 3D boolean array
        connectivity: 6 (face), 18 (face+edge), or 26 (face+edge+corner)
    
    Returns:
        labeled: Array same shape as input with component labels (0=background)
        n_components: Number of components found
    """
    if connectivity == 6:
        structure = ndimage.generate_binary_structure(3, 1)
    elif connectivity == 18:
        structure = ndimage.generate_binary_structure(3, 2)
    elif connectivity == 26:
        structure = ndimage.generate_binary_structure(3, 3)
    else:
        raise ValueError(f"connectivity must be 6, 18, or 26, got {connectivity}")
    
    labeled, n_components = ndimage.label(binary_mask, structure=structure)
    return labeled, n_components


def get_largest_component(labeled: np.ndarray, n_components: int) -> np.ndarray:
    """
    Extract mask of the largest connected component.
    
    Args:
        labeled: Output from connected_components_3d
        n_components: Number of components
    
    Returns:
        Binary mask of largest component only
    """
    if n_components == 0:
        return np.zeros_like(labeled, dtype=bool)
    
    # Count voxels in each component (excluding background=0)
    component_sizes = np.bincount(labeled.ravel())
    component_sizes[0] = 0  # Ignore background
    
    largest_label = np.argmax(component_sizes)
    return labeled == largest_label


def get_component_stats(labeled: np.ndarray, n_components: int) -> list[dict]:
    """
    Get statistics for all connected components.
    
    Returns:
        List of dicts with {label, size, bbox_min, bbox_max}
    """
    stats = []
    
    for label in range(1, n_components + 1):
        mask = labeled == label
        size = int(mask.sum())
        
        if size == 0:
            continue
        
        indices = np.argwhere(mask)
        bbox_min = indices.min(axis=0)
        bbox_max = indices.max(axis=0)
        
        stats.append({
            "label": int(label),
            "size": size,
            "bbox_min": bbox_min.tolist(),
            "bbox_max": bbox_max.tolist(),
        })
    
    # Sort by size descending
    stats.sort(key=lambda x: x["size"], reverse=True)
    return stats


def morphological_clean(
    binary_mask: np.ndarray, 
    operation: str = "opening",
    iterations: int = 1
) -> np.ndarray:
    """
    Apply morphological operations to clean up a binary mask.
    
    Args:
        binary_mask: 3D boolean array
        operation: "opening" (remove small bits), "closing" (fill small holes), 
                  "erosion", "dilation"
        iterations: Number of times to apply operation
    
    Returns:
        Cleaned binary mask
    """
    if operation == "opening":
        return ndimage.binary_opening(binary_mask, iterations=iterations)
    elif operation == "closing":
        return ndimage.binary_closing(binary_mask, iterations=iterations)
    elif operation == "erosion":
        return ndimage.binary_erosion(binary_mask, iterations=iterations)
    elif operation == "dilation":
        return ndimage.binary_dilation(binary_mask, iterations=iterations)
    else:
        raise ValueError(f"Unknown operation: {operation}")