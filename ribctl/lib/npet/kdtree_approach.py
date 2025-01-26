import os
from typing import Tuple
import open3d as o3d
import numpy as np
from data.asset_manager import StructureAssets
from mesh_generation.tunnel_bbox_ptc_constriction import (
    filter_residues_parallel,
    ribosome_entities,
)
from scipy.spatial import cKDTree
import sys

from ribctl.lib.npet.util import landmark_constriction_site, landmark_ptc
data_dir = os.getenv('DATA_DIR')

sys.dont_write_bytecode = True
import pyvista as pv
from npet.mesh_visualization import (
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
    visualize_mesh,
    visualize_pointcloud,
)
from mesh_generation.mesh_full_pipeline import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from mesh_generation.mesh_libsurf import (
    apply_poisson_reconstruction,
    estimate_normals,
    ptcloud_convex_hull_points,
)
# from mesh_generation.util import landmark_constriction_site, landmark_ptc

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

def visualize_clipping_result(
    original_points: np.ndarray,
    clipped_points: np.ndarray,
    mesh_path: str,
    show_mesh: bool = True,
):
    """
    Visualizes the original points, clipped points, and the clipping mesh.

    Parameters:
        original_points (np.ndarray): The original point cloud
        clipped_points (np.ndarray): The clipped point cloud
        mesh_path (str): Path to the mesh file
        show_mesh (bool): Whether to show the mesh or not
    """
    p = pv.Plotter()

    # Add original points in red
    original_cloud = pv.PolyData(original_points)
    p.add_mesh(
        original_cloud,
        color="red",
        point_size=5,
        render_points_as_spheres=True,
        label="Original Points",
    )

    # Add clipped points in blue
    clipped_cloud = pv.PolyData(clipped_points)
    p.add_mesh(
        clipped_cloud,
        color="blue",
        point_size=5,
        render_points_as_spheres=True,
        label="Clipped Points",
    )

    # Add mesh if requested
    if show_mesh:
        mesh = pv.read(mesh_path)
        p.add_mesh(
            mesh, style="wireframe", color="gray", opacity=0.5, label="Clipping Mesh"
        )

    p.add_legend()
    p.show()

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

def create_tunnel_mesh(RCSB_ID:str):
    assets    = StructureAssets(data_dir, RCSB_ID)
    cifpath   = assets.cif_struct
    R         = 30
    H         = 100
    Vsize     = 1
    ATOM_SIZE = 2
    normals_pcd_path       = assets.tunnel_pcd_normal_estimated
    ashape_watertight_mesh = pv.read(assets.ashape_watertight)
    if not os.path.exists(assets.ashape_watertight):
        raise FileNotFoundError(f"File {assets.ashape_watertight} not found")
    # # All atoms
    ptc_pt             = np.array(landmark_ptc(RCSB_ID))
    constriction_pt    = np.array(landmark_constriction_site(RCSB_ID))
    entities           = ribosome_entities(RCSB_ID, cifpath, "R")
    filtered           = filter_residues_parallel(entities, ptc_pt, constriction_pt, R, H)
    filtered_points    = np.array( [atom.get_coord() for residue in filtered for atom in residue.child_list] )
    transformed_points = transform_points_to_C0( filtered_points, ptc_pt, constriction_pt )

    mask, (x, y, z) = create_point_cloud_mask(
        transformed_points,
        radius=R,
        height=H,
        voxel_size=Vsize,
        radius_around_point=ATOM_SIZE,
    )

    points = np.where(~mask)
    empty_coordinates = np.column_stack((x[points[0]], y[points[1]], z[points[2]]))

    back_projected = transform_points_from_C0(
        empty_coordinates, ptc_pt, constriction_pt
    )

    if not os.path.exists(assets.ashape_watertight):
        raise FileNotFoundError(f"File {assets.ashape_watertight} not found")
    select                 = pv.PolyData(back_projected).select_enclosed_points(ashape_watertight_mesh)
    mask                   = select["SelectedPoints"]
    interior               = back_projected[mask == 1]
    empty_in_world_coords = np.array(interior)

    _u_EPSILON_initial_pass     = 5.5
    _u_MIN_SAMPLES_initial_pass = 600
    _u_EPSILON_refinement       = 3.5
    _u_MIN_SAMPLES_refinement   = 175

    d3d_alpha   = 2
    d3d_tol     = 1
    PR_depth    = 6
    PR_ptweight = 3

    db, clusters_container = DBSCAN_capture( empty_in_world_coords, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass )
    #! [ Extract the largest cluster from the DBSCAN clustering ]
    largest_cluster, largest_cluster_id = DBSCAN_pick_largest_cluster(clusters_container)
    # visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs( clusters_container, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass, ptc_pt, constriction_pt, np.array([[1,1,1]]))
    # exit()
    db             , refined_clusters_container = DBSCAN_capture( largest_cluster, _u_EPSILON_refinement, _u_MIN_SAMPLES_refinement )
    refined_cluster, refined_cluster_id         = DBSCAN_pick_largest_cluster(refined_clusters_container)

    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
        clusters_container,
        _u_EPSILON_initial_pass,
        _u_MIN_SAMPLES_initial_pass,
        ptc_pt,
        constriction_pt,
        refined_cluster,
        R,
        H,
    )

    #! [ Transform the cluster back into original coordinate frame ]
    surface_pts = ptcloud_convex_hull_points(refined_cluster, d3d_alpha, d3d_tol)
    visualize_pointcloud(surface_pts, RCSB_ID)

    #! [ Transform the cluster back into Original Coordinate Frame ]
    normal_estimated_pcd = estimate_normals(
        surface_pts,
        kdtree_radius=10,
        kdtree_max_nn=15,
        correction_tangent_planes_n=10,
    )
    o3d.io.write_point_cloud(normals_pcd_path, normal_estimated_pcd)
    apply_poisson_reconstruction(
        normals_pcd_path,
        assets.tunnel_mesh,
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    visualize_mesh(assets.tunnel_mesh)
