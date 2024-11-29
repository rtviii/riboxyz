import numpy as np
from scipy.spatial import cKDTree
import pyvista as pv

from __scripts.cylinder import get_transformation_to_C0, transform_points_to_C0
from mesh_generation.mes_visualization import visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs, visualize_mesh, visualize_pointcloud
from mesh_generation.mesh_full_pipeline import DBSCAN_capture, DBSCAN_pick_largest_cluster
from mesh_generation.mesh_libsurf import apply_poisson_reconstruction, estimate_normals, ptcloud_convex_hull_points
from ribctl.lib.landmarks.constriction import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import filter_residues_parallel, ribosome_entities

from mesh_generation.mesh_paths import convex_hull_cluster_path, surface_with_normals_path, poisson_recon_path

def generate_voxel_centers(radius: float, height: float, voxel_size: float) -> tuple:

    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x  = np.linspace(-radius, radius, nx)
    y  = np.linspace(-radius, radius, ny)
    z  = np.linspace(0, height, nz)
    
    X, Y, Z       = np.meshgrid(x, y, z, indexing='ij')
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)

def create_point_cloud_mask(points: np.ndarray, 
                          radius: float, 
                          height: float,
                          voxel_size: float = 1.0,
                          radius_around_point: float = 2.0):

    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(radius, height, voxel_size)
    tree = cKDTree(points)
    indices = tree.query_ball_point(voxel_centers, radius_around_point)
    

    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True
    
    point_cloud_mask = point_cloud_mask.reshape(grid_shape)
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    cylinder_mask = (np.sqrt(X**2 + Y**2) <= radius)
    hollow_cylinder = ~cylinder_mask
    
    final_mask = hollow_cylinder | point_cloud_mask
    return final_mask, (x, y, z)

def transform_points_to_C0(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_translated  = points + translation
    points_transformed = points_translated @ rotation.T
    
    return points_transformed

def transform_points_from_C0(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray) -> np.ndarray:

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

def verify_mesh_quality(mesh_path: str) -> dict:
    """
    Verifies the quality of the input mesh and returns diagnostics.
    """
    mesh = pv.read(mesh_path)
    stats = {
        'n_points': mesh.n_points,
        'n_faces': mesh.n_faces,
        'is_manifold': mesh.is_manifold,
        'bounds': mesh.bounds,
    }
    
    try:
        stats['volume'] = mesh.volume
    except:
        stats['volume'] = None
        print("Warning: Could not compute mesh volume")
        
    return stats

def visualize_clipping_result(original_points: np.ndarray, 
                            clipped_points: np.ndarray, 
                            mesh_path: str,
                            show_mesh: bool = True):
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
    p.add_mesh(original_cloud, color='red', point_size=5, render_points_as_spheres=True,
              label='Original Points')
    
    # Add clipped points in blue
    clipped_cloud = pv.PolyData(clipped_points)
    p.add_mesh(clipped_cloud, color='blue', point_size=5, render_points_as_spheres=True,
              label='Clipped Points')
    
    # Add mesh if requested
    if show_mesh:
        mesh = pv.read(mesh_path)
        p.add_mesh(mesh, style='wireframe', color='gray', opacity=0.5, label='Clipping Mesh')
    
    p.add_legend()
    p.show()



def main():
    # Load your points and transform them as before
    RCSB_ID    = '4ug0'.upper()
    R          = 40
    H          = 100
    Vsize      = 1
    ATOM_SIZE  = 2
    base_point = np.array(PTC_location(RCSB_ID).location)
    axis_point = np.array(get_constriction(RCSB_ID) )

    residues           = filter_residues_parallel( ribosome_entities(RCSB_ID, 'R'), base_point, axis_point, R, H)
    points             = np.array([atom.get_coord() for residue in residues for atom in residue.child_list])
    np.save('LSU_points.npy', points)

    transformed_points = transform_points_to_C0(points, base_point, axis_point)

    mask, (x, y, z) = create_point_cloud_mask(
        transformed_points,
        radius              = R,
        height              = H,
        voxel_size          = Vsize,
        radius_around_point = ATOM_SIZE
    )
    
    points = np.where(~mask)
    empty_coordinates = np.column_stack((
        x[points[0]], 
        y[points[1]], 
        z[points[2]]
    ))
    back_projected    = transform_points_from_C0(empty_coordinates ,base_point,axis_point)

    # occupied_points = pv.PolyData(empty_coordinates)
    # world_coords_pv = pv.PolyData(empty_in_world_coords)
    # visualize_pointcloud(occupied_points, world_coords_pv)

    mesh_path = "alpha_shape_watertight.ply"
    mesh      = pv.read(mesh_path)
    select    = pv.PolyData( back_projected ).select_enclosed_points(mesh)
    mask      = select['SelectedPoints']
    interior  = back_projected[mask == 1]

    empty_in_world_coords = np.array(interior)



    _u_EPSILON_initial_pass = 5.5 
    _u_MIN_SAMPLES_initial_pass =  600
    _u_EPSILON_refinement = 3.5
    _u_MIN_SAMPLES_refinement = 175

    d3d_alpha = 2
    d3d_tol = 1
    PR_depth = 6
    PR_ptweight = 3
    db, clusters_container = DBSCAN_capture(empty_in_world_coords , _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass)
    #! [ Extract the largest cluster from the DBSCAN clustering ]
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs( clusters_container, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass)
    largest_cluster = DBSCAN_pick_largest_cluster(clusters_container)
    db, clusters_refinement = DBSCAN_capture(largest_cluster , _u_EPSILON_refinement, _u_MIN_SAMPLES_refinement)
    refined = DBSCAN_pick_largest_cluster(clusters_refinement)


    #! [ Transform the cluster back into original coordinate frame ]
    surface_pts = ptcloud_convex_hull_points( refined, d3d_alpha, d3d_tol )
    visualize_pointcloud( surface_pts, RCSB_ID, False, "{}.surface_pts.gif".format(RCSB_ID) )

    #! [ Transform the cluster back into Original Coordinate Frame ]
    np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    estimate_normals(
        surface_pts,
        surface_with_normals_path(RCSB_ID),
        kdtree_radius               = 10,
        kdtree_max_nn               = 15,
        correction_tangent_planes_n = 10,
    )
    apply_poisson_reconstruction(
        surface_with_normals_path(RCSB_ID),
        poisson_recon_path(RCSB_ID),
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    visualize_mesh(
        pv.read(poisson_recon_path(RCSB_ID)),
        RCSB_ID
    )

    

if __name__ == '__main__':
    main()


if __name__ == "__main__":
    # Example points (replace with your actual points)
    points = np.load('tunnel_points.npy')
    print(points.shape)
    
    mesh_path      = "alpha_shape_watertight.ply"
    mesh           = pv.read(mesh_path)
    points_poly = pv.PolyData(points)
    select = points_poly.select_enclosed_points(mesh)
    print(points.shape)
    mask = select['SelectedPoints']
    interior = points[mask == 1]
    exterior = points[mask == 0]
    pl = pv.Plotter()
    _ = pl.add_mesh(mesh, style='wireframe')
    _ = pl.add_points(interior, color='r')
    _ = pl.add_points(exterior, color='b')
    pl.show()