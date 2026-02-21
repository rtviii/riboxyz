# npet2/backends/geometry.py
"""
Geometric primitives ported from riboxyz's kdtree_approach.py and alphalib.py.

Only the functions that npet2 actually uses are included here.
"""
from __future__ import annotations

import multiprocessing
import subprocess
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Callable, List, Optional, Tuple, TypeVar

import numpy as np
import open3d as o3d
import plyfile
import pyvista as pv
from Bio.PDB.MMCIFParser import MMCIFParser
from scipy.spatial import cKDTree
from sklearn.cluster import DBSCAN as skDBSCAN

from npet2.core.config import SETTINGS

T = TypeVar("T")


# ---------------------------------------------------------------------------
# Coordinate transforms (PTC/constriction axis -> canonical C0 frame)
# ---------------------------------------------------------------------------

def get_transformation_to_C0(
    base_point: np.ndarray, axis_point: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    z_axis = np.array([0, 0, 1])

    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])
    else:
        v = np.cross(axis_unit, z_axis)
        s = np.linalg.norm(v)
        c = np.dot(axis_unit, z_axis)
        v_skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + v_skew + (v_skew @ v_skew) * (1 - c) / (s * s)

    return -base_point, R


def transform_points_to_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    return (points + translation) @ rotation.T


def transform_points_from_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    return (points @ rotation) - translation


# ---------------------------------------------------------------------------
# DBSCAN wrappers
# ---------------------------------------------------------------------------

def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps: float,
    min_samples: int,
    metric: str = "euclidean",
) -> Tuple[skDBSCAN, dict]:
    print(
        f"Running DBSCAN on {len(ptcloud)} points. "
        f"eps={eps}, min_samples={min_samples}, distance_metric={metric}"
    )
    db = skDBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(ptcloud)

    clusters: dict[int, list] = {}
    for point, label in zip(ptcloud, db.labels_):
        clusters.setdefault(int(label), []).append(point)

    return db, dict(sorted(clusters.items()))


def DBSCAN_pick_largest_cluster(
    clusters_container: dict[int, list],
) -> Tuple[np.ndarray, int]:
    best_id = 0
    for k, v in clusters_container.items():
        if int(k) == -1:
            continue
        if len(v) > len(clusters_container.get(best_id, [])):
            best_id = int(k)
    return np.array(clusters_container[best_id]), best_id


# ---------------------------------------------------------------------------
# Poisson reconstruction (calls external binary)
# ---------------------------------------------------------------------------

def apply_poisson_reconstruction(
    surf_estimated_ptcloud_path: str,
    output_path: Path,
    recon_depth: int = 6,
    recon_pt_weight: float = 3.0,
) -> None:
    print(f"Rolling Poisson Reconstruction: {surf_estimated_ptcloud_path} -> {output_path}")
    command = [
        SETTINGS.poisson_recon_bin,
        "--in", str(surf_estimated_ptcloud_path),
        "--out", str(output_path),
        "--depth", str(recon_depth),
        "--pointWeight", str(recon_pt_weight),
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        data = plyfile.PlyData.read(str(output_path))
        data.text = True
        ascii_path = str(output_path).rsplit(".", 1)[0] + "_ascii.ply"
        data.write(ascii_path)
        print(f">>Wrote {output_path} and the _ascii version.")
    else:
        print(f">>Error: {process.stderr}")


# ---------------------------------------------------------------------------
# Cylinder filtering
# ---------------------------------------------------------------------------

def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    z_min: float = 0.0,
) -> bool:
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    point_vector = point - base_point
    projection = np.dot(point_vector, axis_unit)

    projection_point = base_point + projection * axis_unit
    radial_distance = np.linalg.norm(point - projection_point)

    return (radial_distance <= radius) and (z_min <= projection <= height)


def _worker_process_chunk(chunk_data):
    positions, base_point, axis_point, radius, height, indices, z_min = chunk_data
    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(pos, base_point, axis_point, radius, height, z_min=z_min):
            results.append(indices[i])
    return results


def get_residue_position(residue):
    return residue.center_of_mass()


def filter_residues_parallel(
    residues: list,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None,
    z_min: float = 0.0,
) -> list:
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    if chunk_size is None:
        chunk_size = max(1, len(residues) // (max_workers * 4))

    positions = np.array([get_residue_position(r) for r in residues])
    indices = list(range(len(residues)))
    index_chunks = [indices[i:i + chunk_size] for i in range(0, len(indices), chunk_size)]

    chunks_data = [
        (positions[idx], base_point, axis_point, radius, height, idx, z_min)
        for idx in index_chunks
    ]

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))

    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]


# ---------------------------------------------------------------------------
# Grid voxel mask (legacy kdtree backend)
# ---------------------------------------------------------------------------

def generate_voxel_centers(
    radius: float, height: float, voxel_size: float, z_min: float = 0.0
) -> Tuple[np.ndarray, tuple]:
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(z_min, z_min + height, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)


def create_point_cloud_mask(
    points: np.ndarray,
    radius: float,
    height: float,
    voxel_size: float = 1.0,
    radius_around_point: float = 2.0,
    z_min: float = 0.0,
):
    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(
        radius, height, voxel_size, z_min=z_min
    )
    tree = cKDTree(points)
    indices = tree.query_ball_point(voxel_centers, radius_around_point)

    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True
    point_cloud_mask = point_cloud_mask.reshape(grid_shape)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    cylinder_mask = np.sqrt(X**2 + Y**2) <= radius
    hollow_cylinder = ~cylinder_mask

    final_mask = hollow_cylinder | point_cloud_mask
    return final_mask, (x, y, z)


# ---------------------------------------------------------------------------
# Normal estimation
# ---------------------------------------------------------------------------

def estimate_normals(
    convex_hull_surface_pts: np.ndarray,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
) -> o3d.geometry.PointCloud:
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    return pcd


# ---------------------------------------------------------------------------
# Alpha shape / surface extraction (from alphalib)
# ---------------------------------------------------------------------------

def cif_to_point_cloud(
    cif_path: str,
    chains: list[str] | None = None,
    do_atoms: bool = False,
    exclude_chains: list[str] = [],
) -> np.ndarray:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)
    coordinates = []

    first_model = structure[0]
    if do_atoms:
        for chain in first_model:
            if chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())
    else:
        for chain in first_model:
            if chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                coordinates.append(residue.center_of_mass())

    if not coordinates:
        raise ValueError(f"No coordinates found in {cif_path}")

    return np.array(coordinates)


def quick_surface_points(
    pointcloud: np.ndarray, alpha: float, tolerance: float, offset: float
) -> np.ndarray:
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=alpha, tol=tolerance, offset=offset, progress_bar=True)
    surface = grid.extract_surface().cast_to_pointset()
    return surface.points


def fast_normal_estimation(
    surface_pts: np.ndarray,
    kdtree_radius: float,
    max_nn: int,
    tangent_planes_k: int,
) -> o3d.geometry.PointCloud:
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(surface_pts)
    search_param = o3d.geometry.KDTreeSearchParamHybrid(radius=kdtree_radius, max_nn=max_nn)
    pcd.estimate_normals(search_param=search_param)
    pcd.orient_normals_consistent_tangent_plane(k=tangent_planes_k)
    return pcd


def validate_mesh_pyvista(mesh_or_path, stage="unknown") -> bool:
    if isinstance(mesh_or_path, (str, Path)):
        import os
        if not os.path.exists(mesh_or_path):
            print(f"WARNING: Mesh file does not exist: {mesh_or_path}")
            return False
        try:
            mesh = pv.read(str(mesh_or_path))
        except Exception as e:
            print(f"WARNING: Failed to load mesh at {mesh_or_path}: {e}")
            return False
    else:
        mesh = mesh_or_path

    if mesh is None:
        print(f"WARNING: Null mesh at stage {stage}")
        return False

    try:
        edges = mesh.extract_feature_edges(
            boundary_edges=True,
            feature_edges=False,
            manifold_edges=False,
            non_manifold_edges=False,
        )
        return edges.n_cells == 0
    except Exception as e:
        print(f"WARNING: Error checking watertightness: {e}")
        return False