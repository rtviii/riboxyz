# ribctl/lib/npet2/backends/meshing.py
from __future__ import annotations
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.ndimage import binary_fill_holes, gaussian_filter

def save_mesh_with_ascii(mesh: pv.PolyData, path: Path, tag: str = "") -> None:
    """Save mesh as binary PLY + ASCII PLY side by side."""
    mesh.save(str(path))
    ascii_path = path.parent / f"{path.stem}_ascii.ply"
    try:
        mesh.save(str(ascii_path), binary=False)
    except Exception:
        try:
            import plyfile
            data = plyfile.PlyData.read(str(path))
            data.text = True
            data.write(str(ascii_path))
        except Exception as e:
            print(f"[meshing] ASCII PLY write failed{' (' + tag + ')' if tag else ''}: {e}")

def mesh_from_binary_volume(
    mask: np.ndarray,
    origin: np.ndarray,
    voxel_size: float,
    *,
    gaussian_sigma_voxels: float = 1.5,
    smooth_method: str = "taubin",
    smooth_iters: int = 20,
    taubin_pass_band: float = 0.1,
    fill_holes_size: float = 100.0,
) -> tuple[pv.PolyData, pv.PolyData]:
    """
    Marching cubes on a Gaussian-blurred binary volume, followed by mesh smoothing.

    Returns (smoothed_mesh, pre_smooth_mesh).
    Both are in the same coordinate frame as the input volume.
    Caller is responsible for coordinate transforms and saving.
    """
    vol = np.pad(mask.astype(np.float32), 2, constant_values=0.0)
    origin_pad = np.asarray(origin, dtype=np.float32) - 2 * voxel_size

    if gaussian_sigma_voxels > 0:
        vol = gaussian_filter(vol, sigma=gaussian_sigma_voxels)

    img = pv.ImageData(
        dimensions=vol.shape,
        spacing=(voxel_size, voxel_size, voxel_size),
        origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
    )
    img.point_data["values"] = vol.ravel(order="F")

    surf = img.contour(isosurfaces=[0.5], scalars="values").triangulate()
    if surf.n_points == 0:
        raise ValueError("Marching cubes produced empty surface")

    surf = surf.clean(tolerance=0.0)

    if fill_holes_size > 0:
        surf = surf.fill_holes(fill_holes_size)

    surf = surf.connectivity(largest=True)

    pre_smooth = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)

    if smooth_iters > 0:
        if smooth_method == "taubin":
            surf = surf.smooth_taubin(n_iter=smooth_iters, pass_band=taubin_pass_band)
        else:
            surf = surf.smooth(n_iter=smooth_iters)

    surf = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)
    return surf, pre_smooth


def voxelize_points(
    points: np.ndarray,
    voxel_size: float,
    pad_voxels: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    lo = points.min(axis=0) - pad_voxels * voxel_size
    hi = points.max(axis=0) + pad_voxels * voxel_size

    shape = tuple(np.ceil((hi - lo) / voxel_size).astype(int) + 1)

    ijk = np.floor((points - lo) / voxel_size + 0.5).astype(np.int32)
    for d in range(3):
        ijk[:, d] = np.clip(ijk[:, d], 0, shape[d] - 1)

    mask = np.zeros(shape, dtype=bool)
    mask[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True

    mask = binary_fill_holes(mask)
    return mask, lo.astype(np.float32)


def clip_mesh_to_atom_clearance(
    mesh: pv.PolyData,
    atom_xyz: np.ndarray,
    min_clearance_A: float = 1.5,
) -> pv.PolyData:
    from scipy.spatial import cKDTree

    tree = cKDTree(atom_xyz)
    pts = np.asarray(mesh.points, dtype=np.float64)

    dist, idx = tree.query(pts, k=1)
    violating = dist < min_clearance_A

    if violating.sum() == 0:
        return mesh

    nearest = atom_xyz[idx[violating]]
    direction = pts[violating] - nearest
    norms = np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-8)
    pts[violating] = nearest + (direction / norms) * min_clearance_A

    result = mesh.copy()
    result.points = pts.astype(np.float32)
    return result