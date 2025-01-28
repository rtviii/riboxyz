from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os
from typing import Callable, List, Literal, Optional, Tuple, TypeVar
import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import sys

data_dir = os.getenv('DATA_DIR')
sys.dont_write_bytecode = True
from Bio.PDB.Entity import Entity
from Bio.PDB.MMCIFParser import FastMMCIFParser
import numpy as np
from Bio.PDB.Chain import Chain
import pyvista as pv
from pathlib import Path
from pprint import pprint
import subprocess
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import numpy as np
import plyfile
import warnings
from ribctl import POISSON_RECON_BIN
import numpy as np
from sklearn.cluster import DBSCAN
import requests
warnings.filterwarnings("ignore")
import os

#
# POISSON_RECON_BIN = os.getenv("POISSON_RECON_BIN")


def landmark_constriction_site(rcsb_id: str) -> np.ndarray:
    """
    Fetches the constriction site location for a given RCSB PDB ID.
    
    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")
        
    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the constriction site
    """
    url = f"http://localhost:8000/structures/constriction_site?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}
    
    response = requests.get(url, headers=headers)
    response.raise_for_status()  # Raises an exception for 4XX/5XX status codes
    
    data = response.json()
    return np.array(data["location"])

def landmark_ptc(rcsb_id: str) -> np.ndarray:
    """
    Fetches the peptidyl transferase center (PTC) location for a given RCSB PDB ID.
    
    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")
        
    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the PTC site
    """
    url = f"http://localhost:8000/structures/ptc?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}
    
    response = requests.get(url, headers=headers)
    response.raise_for_status()
    
    data = response.json()
    return np.array(data["location"])


def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps           ,
    min_samples   ,
    metric        : str = "euclidean",
): 

    u_EPSILON     = eps
    u_MIN_SAMPLES = min_samples
    u_METRIC      = metric

    print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(ptcloud), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) ) 
    db     = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit( ptcloud )
    labels = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(ptcloud, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))
    return db, CLUSTERS_CONTAINER

def DBSCAN_pick_largest_cluster(clusters_container:dict[int,list], pick_manually:bool=False)->tuple[np.ndarray, int]:
    DBSCAN_CLUSTER_ID = 0
    # if pick_manually:
    #     print("-------------------------------")
    #     print("Running Manual Cluster Selection")
    #     picked_id =  int(input("Enter Cluster ID to proceed the reconstruction with\n (options:[{}]):".format(list( clusters_container.keys() ))))
    #     print("Choise cluster # {}".format(picked_id))
    #     if picked_id == -2:
    #         # if picked -2 ==> return largest
    #         for k, v in clusters_container.items():
    #             if int(k) == -1:
    #                 continue
    #             elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
    #                 DBSCAN_CLUSTER_ID = int(k)
    #         return np.array(clusters_container[DBSCAN_CLUSTER_ID])

    #     return np.array(clusters_container[picked_id])

    for k, v in clusters_container.items():
        # print(f"Cluster {k} has {len(v)} points.")
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)
    return np.array(clusters_container[DBSCAN_CLUSTER_ID]), DBSCAN_CLUSTER_ID


def apply_poisson_reconstruction(
    surf_estimated_ptcloud_path: str,
    output_path: Path,
    recon_depth: int = 6,
    recon_pt_weight: int = 3,
):
    # The documentation can be found at https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version16.04/ in "PoissonRecon" binary
    command = [
        POISSON_RECON_BIN,
        "--in",
        surf_estimated_ptcloud_path,
        "--out",
        output_path,
        "--depth",
        str(recon_depth),
        "--pointWeight",
        str(recon_pt_weight),
        "--threads 8",
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        data = plyfile.PlyData.read(output_path)
        data.text = True
        ascii_duplicate = output_path.as_posix().split(".")[0] + "_ascii.ply"
        data.write(ascii_duplicate)
        print(">>Wrote {} and the _ascii version.".format(output_path))
    else:
        print(">>Error:", process.stderr)

def ptcloud_convex_hull_points(
    pointcloud: np.ndarray, ALPHA: float, TOLERANCE: float
) -> np.ndarray:
    assert pointcloud is not None
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points

def estimate_normals(
    convex_hull_surface_pts: np.ndarray,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    # o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    return pcd

T = TypeVar('T')

def generate_voxel_centers(radius: float, height: float, voxel_size: float) -> tuple:
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)

def create_point_cloud_mask(
    points: np.ndarray,
    radius: float,
    height: float,
    voxel_size: float = 1.0,
    radius_around_point: float = 2.0,
):

    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(
        radius, height, voxel_size
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

def get_transformation_to_C0(
    base_point: np.ndarray, axis_point: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute transformation matrices to move arbitrary cylinder to C0 configuration.
    Returns translation vector and rotation matrix.
    """
    # Get cylinder axis vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    # Get rotation that aligns axis_unit with [0, 0, 1]
    z_axis = np.array([0, 0, 1])

    # Use Rodrigues rotation formula to find rotation matrix
    # that rotates axis_unit to z_axis
    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])  # 180-degree rotation around x-axis
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
    points_translated = points + translation
    points_transformed = points_translated @ rotation.T

    return points_transformed

def transform_points_from_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:

    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_unrotated = points @ rotation
    points_untranslated = points_unrotated - translation

    return points_untranslated

def clip_pointcloud_with_mesh(points: np.ndarray, mesh_path: str) -> np.ndarray:
    """
    Clips a point cloud to keep only points that lie inside a mesh using point-by-point checking.

    Parameters:
        points (np.ndarray): Nx3 array of points to clip
        mesh_path (str): Path to the mesh file (PLY format)

    Returns:
        np.ndarray: Points that lie inside the mesh
    """
    # Load the mesh
    mesh = pv.read(mesh_path)

    # Convert to PolyData if needed
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()

    # Ensure mesh is triangulated
    if not mesh.is_all_triangles:
        mesh = mesh.triangulate()
    else:
        print("Mesh is triangulated")

    # Initialize mask array
    mask = np.zeros(len(points), dtype=bool)
    print("Got mask", mask.shape)

    # Check each point
    for i, point in enumerate(points):
        # Create a single-point PolyData
        point_data = pv.PolyData(point.reshape(1, 3))
        # Check if point is inside
        selection = mesh.select_enclosed_points(point_data, check_surface=False)
        mask[i] = bool(selection.point_data["SelectedPoints"][0])

        # Progress indicator
        if i % 1000 == 0:
            print(f"Processed {i}/{len(points)} points")

    # Return filtered points
    clipped_points = points[mask]
    print(f"Kept {len(clipped_points)}/{len(points)} points")

    return clipped_points

def verify_mesh_quality(mesh) -> dict:
    """
    Verifies the quality of the input mesh and returns diagnostics.
    """
    stats = {
        "n_points": mesh.n_points,
        "n_faces": mesh.n_faces,
        "is_manifold": mesh.is_manifold,
        "bounds": mesh.bounds,
        "open_edges": mesh.n_open_edges,
    }

    try:
        stats["volume"] = mesh.volume
    except:
        stats["volume"] = None
        print("Warning: Could not compute mesh volume")

    return stats

def clip_pcd_via_ashape(
    pcd: np.ndarray, mesh: pv.PolyData
) -> Tuple[np.ndarray, np.ndarray]:

    # TODO
    points_poly = pv.PolyData(pcd)
    select = points_poly.select_enclosed_points(mesh)
    mask = select["SelectedPoints"]
    ashape_interior = pcd[mask == 1]
    ashape_exterior = pcd[mask == 0]
    return ashape_interior, ashape_exterior

def ribosome_entities(rcsb_id:str, cifpath:str, level=Literal['R']|Literal[ 'A' ], skip_nascent_chain:List[str]=[])->list[Entity]:
    structure = FastMMCIFParser(QUIET=True).get_structure(rcsb_id, cifpath)
    residues = []
    for chain in structure.child_list[0]:
        if  chain.id in skip_nascent_chain:
            print("Skipping nascent chain ", chain.id)
            continue
        residues.extend(chain)
    # [residues.extend(chain) for chain in structure.child_list[0] ]
    if level == 'R':
        return residues
    elif level == 'A':
        atoms = []
        [atoms.extend(a) for a in residues]
        return atoms
    else:
        raise

def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float

) -> bool:

    point      = np.asarray(point)
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)
    
    # Calculate cylinder axis direction vector
    axis        = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit   = axis / axis_length
    
    # Calculate vector from base to point
    point_vector = point - base_point
    
    # Project vector onto cylinder axis
    projection = np.dot(point_vector, axis_unit)
    
    # Calculate perpendicular vector from axis to point
    projection_point = base_point + projection * axis_unit
    radial_vector    = point - projection_point
    
    # Calculate radial distance
    radial_distance = np.linalg.norm(radial_vector)
    
    # Check if point is inside cylinder
    return (radial_distance <= radius) and (0 <= projection <= height)

def make_cylinder_predicate(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray] 
) -> Callable[[T], bool]:
    def predicate(obj: T) -> bool:
        position = position_getter(obj)
        return is_point_in_cylinder(position, base_point, axis_point, radius, height)
    return predicate

def get_residue_position(residue):
    return residue.center_of_mass()

def _worker_process_chunk(chunk_data):
    """
    Worker function that processes a chunk of residues.
    Takes a tuple of (residue_positions, base_point, axis_point, radius, height, indices)
    Returns indices of residues that are inside the cylinder.
    """
    positions, base_point, axis_point, radius, height, indices = chunk_data
    
    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(pos, base_point, axis_point, radius, height):
            results.append(indices[i])
    
    return results

def filter_residues_parallel(
    residues: List[T],
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None
) -> List[T]:
    """
    Filter residues in parallel using ProcessPoolExecutor.
    
    Parameters:
    -----------
    residues : List[T]
        List of residue objects to filter
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    chunk_size : Optional[int]
        Size of chunks to process in parallel. If None, calculated automatically
    max_workers : Optional[int]
        Maximum number of worker processes. If None, uses CPU count
    
    Returns:
    --------
    List[T]
        Filtered list of residues whose positions lie within the cylinder
    """
    # Set defaults for parallel processing parameters
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    if chunk_size is None:
        # Aim for ~4 chunks per worker for better load balancing
        chunk_size = max(1, len(residues) // (max_workers * 4))
    
    # Pre-compute all positions and create index mapping
    positions = np.array([get_residue_position(r) for r in residues])
    
    indices = list(range(len(residues)))
    index_chunks = [
        indices[i:i + chunk_size] 
        for i in range(0, len(indices), chunk_size)
    ]
    
    # Create data chunks for processing
    chunks_data = [ (positions[idx], base_point, axis_point, radius, height, idx) for idx in index_chunks ]
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))
    
    # Flatten results and get corresponding residues
    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]
