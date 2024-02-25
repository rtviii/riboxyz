from pprint import pprint
import typing
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import json
import os
import numpy as np
import numpy as np
from sklearn.cluster import DBSCAN
from mesh_generation.bbox_extraction import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_bbox,
    parse_struct_via_centerline,
)

from compas.geometry import bounding_box
from mesh_generation.voxelize import expand_atomcenters_to_spheres_threadpool, normalize_atom_coordinates
from ribctl import  EXIT_TUNNEL_WORK, RIBETL_DATA

#? ---------- Params ------------
RCSB_ID = "6Z6K"

u_EPSILON     = 5.5
u_MIN_SAMPLES = 600
u_METRIC      = "euclidean"

DBSCAN_CLUSTER_ID = 3
#? ---------- Paths ------------
tunnel_atom_encoding_path    = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_tunnel_atoms_bbox.json".format(rcsb_id))
normalization_vectors_path   = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_normalization_vectors.npy".format(rcsb_id))
selected_dbscan_cluster_path = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_dbscan_cluster.npy".format(rcsb_id))
convex_hull_cluster_path     = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_convex_hull.npy".format(rcsb_id))
surface_with_normals_path    = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_normal_estimated_surf.ply".format(rcsb_id))
poisson_recon_path           = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK, "{}_poisson_recon.ply".format(rcsb_id))
ptc_data_path                = lambda rcsb_id : os.path.join(RIBETL_DATA,rcsb_id, "{}_PTC_COORDINATES.json".format(rcsb_id))
#? ------------------------------


def bounding_box(points:np.ndarray):
    """Computes the axis-aligned minimum bounding box of a list of points.

    Parameters
    ----------
    points : sequence[point]
        XYZ coordinates of the points.

    Returns
    -------
    list[[float, float, float]]
        XYZ coordinates of 8 points defining a box.


    """
    x, y, z = zip(*points)
    min_x = min(x)
    max_x = max(x)
    min_y = min(y)
    max_y = max(y)
    min_z = min(z)
    max_z = max(z)
    return [
        [min_x, min_y, min_z],
        [max_x, min_y, min_z],
        [max_x, max_y, min_z],
        [min_x, max_y, min_z],
        [min_x, min_y, max_z],
        [max_x, min_y, max_z],
        [max_x, max_y, max_z],
        [min_x, max_y, max_z],
    ]

def extract_bbox_atoms(rcsb_id:str)->list:

    centerline_expansion_atoms       = parse_struct_via_centerline(rcsb_id, open_tunnel_csv(rcsb_id))
    centerline_expansion_coordinates = np.array([a.get_coord() for a in centerline_expansion_atoms])
    bbox                             = bounding_box(centerline_expansion_coordinates)
    bbox_interior_atoms              = parse_struct_via_bbox(rcsb_id, bbox)

    return encode_atoms(rcsb_id, bbox_interior_atoms, write=True, writepath=tunnel_atom_encoding_path(RCSB_ID))

def expand_bbox_atoms_to_spheres(bbox_data:list|None ):

    if bbox_data is None:
        with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile:
            bbox_data = json.load(infile)

    __cords = np.array(list(map(lambda x: x["coord"], bbox_data)))
    __radii = np.array(list(map(lambda x: x["vdw_radius"], bbox_data)))

    sphere_sources = zip(__cords, __radii)

    SINK     = []
    expanded = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)


    return np.array(expanded)

def index_grid(expanded_sphere_voxels:np.ndarray):
    voxel_size =1
    normalized_sphere_cords, mean_abs_vectors = normalize_atom_coordinates(expanded_sphere_voxels)
    sphere_cords_quantized  = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)

    max_values      = np.max(sphere_cords_quantized, axis=0)
    grid_dimensions = max_values + 1
    vox_grid        = np.zeros(grid_dimensions)

    vox_grid[
        sphere_cords_quantized[:, 0],
        sphere_cords_quantized[:, 1],
        sphere_cords_quantized[:, 2],
    ] = 1  

    __xyz_v_negative_ix = np.asarray( np.where(vox_grid != 1) )  
    __xyz_v_positive_ix = np.asarray( np.where(vox_grid == 1) )  # get back indexes of populated voxels

    return  __xyz_v_positive_ix.T, __xyz_v_negative_ix.T, grid_dimensions, mean_abs_vectors

def interior_capture_DBSCAN(xyz_v_negative:np.ndarray):

    attempts: list[typing.Tuple[float, int, str]] = [
        ( 2  , 60 , "euclidean"    ),
        ( 4  , 100, "euclidean"  , ), # In 0.84s. Estimated number of [ clusters, noise points ]: [ 3, 2610 ]
        ( 5  , 200, "euclidean"  , ), # In 0.84s. Estimated number of [ clusters, noise points ]: [ 3, 2610 ] (
        ( 6  , 400, "euclidean"    ), # In 1.583s. Estimated number of [ clusters, noise points ]: [ 2, 8442 ]
        ( 5  , 500, "euclidean"  , ), # <<<< Really good result In 0.98s. Estimated number of [ clusters, noise points ]: [ 10, 123078 ]
        ( 6  , 600, "euclidean"  , ), # <<<< Best so far In 1.46s. Estimated number of [ clusters, noise points ]: [ 9, 53038 ]
        ( 5.5, 600, "euclidean"    ), # <<<< A little finer res
        ( 5.5, 650, "euclidean"    ), #
        ( 5.5, 650, "canberra"     ), # failed
        ( 5.5, 650, "sqeuclidean"  ), # 0 clusters
        ( 4  , 100, "sqeuclidean"  ),
        ( 3  , 40 , "canberra"     ),
        ( 6  , 600, "euclidean"  , ), # <<<< Best so far In 1.46s. Estimated number of [ clusters, noise points ]: [ 9, 53038 ]
        ( 5.5, 600, "euclidean"    ), # <<<< A little finer res
    ]

    cluster_colors = dict(zip(range(-1, 40), plt.cm.terrain(np.linspace(0, 1, 40))))
    for k, v in cluster_colors.items():
        cluster_colors[k] = [*v[:3], 0.5]


    u_EPSILON     = attempts[-1][0]
    u_MIN_SAMPLES = attempts[-1][1]
    u_METRIC      = attempts[-1][2]

    print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(xyz_v_negative), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) )
    db                 = DBSCAN(eps=u_EPSILON, min_samples=u_MIN_SAMPLES, metric=u_METRIC, n_jobs=5).fit( xyz_v_negative )
    labels             = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(xyz_v_negative, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))


    return db, CLUSTERS_CONTAINER

def DBSCAN_CLUSTERS_visualize_all(dbscan_cluster_dict:dict[int, list]):

    for k, v in dbscan_cluster_dict.items():
        print("Cluster {} has {} points.".format(k, len(v)))

    clusters_palette = dict(zip(range(-1, 40), plt.cm.terrain(np.linspace(0, 1, 40))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]

    combined_cluster_colors = []
    combined_cluster_points = []
    for dbscan_label, coordinates in dbscan_cluster_dict.items():
        combined_cluster_points.extend(coordinates)
        combined_cluster_colors.extend([ clusters_palette[dbscan_label * 2] if dbscan_label != -1 else [0, 0, 0, 1]]*len(coordinates))

    point_cloud         = pv.PolyData(combined_cluster_points)
    point_cloud["rgba"] = combined_cluster_colors
    point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)

def DBSCAN_CLUSTERS_visualize_one(positive_space:np.ndarray, selected_cluster:np.ndarray):
    rgbas_cluster  = [[15, 10, 221, 1] for datapoint in selected_cluster]
    rgbas_positive = np.array([[205, 209, 228, 0.2] for _ in positive_space])
    combined       = np.concatenate([selected_cluster, positive_space])
    rgbas_combined = np.concatenate([rgbas_cluster, rgbas_positive])

    point_cloud         = pv.PolyData(combined)
    point_cloud["rgba"] = rgbas_combined
    point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)

def surface_pts_via_convex_hull(selected_cluster:np.ndarray|None=None):
    if not os.path.exists(convex_hull_cluster_path(RCSB_ID)):
        assert(selected_cluster is not None)
        cloud       = pv.PolyData(selected_cluster)
        grid        = cloud.delaunay_3d(alpha=2, tol=1.5,offset=2,progress_bar=True)
        convex_hull = grid.extract_surface().cast_to_pointset()
        np.save(convex_hull_cluster_path(RCSB_ID), convex_hull.points)
        print("Saved generated convex hull to {}".format(convex_hull_cluster_path(RCSB_ID)))
        return convex_hull.points
    else:
        pts = np.load(convex_hull_cluster_path(RCSB_ID))
        print("Loaded generated convex hull from {}".format(convex_hull_cluster_path(RCSB_ID)))
        return pts
    # convex_hull.plot(show_edges=True)

def estimate_normals(convex_hull_surface_pts:np.ndarray|None=None):
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=15))
    pcd.orient_normals_consistent_tangent_plane(k=10)
    o3d.visualization.draw_geometries([pcd], point_show_normal=True, window_name="Point Cloud with Normals")
    o3d.io.write_point_cloud(surface_with_normals_path(RCSB_ID), pcd)
    print("Wrote surface with normals {}".format(surface_with_normals_path(RCSB_ID)))


def move_cords_to_normalized_cord_frame(grid_dimensions:np.ndarray, translation_vectors:np.ndarray, original_cords:np.ndarray):
    """this is a helper function for plotting to move additional atom coordinates into the cord frame of the mesh (vox grid indices)"""
    normalized_original_cords           = original_cords - translation_vectors[0] + translation_vectors[1]
    voxel_size                          = 1
    normalized_original_cords_quantized = np.round( normalized_original_cords / voxel_size).astype(int)
    vox_grid                            = np.zeros(grid_dimensions)
    vox_grid[
        normalized_original_cords_quantized[:, 0],
        normalized_original_cords_quantized[:, 1],
        normalized_original_cords_quantized[:, 2]
        ] = 1 
    __xyz_v_positive_ix = np.asarray( np.where(vox_grid == 1) )  # get back indexes of populated voxels
    return  __xyz_v_positive_ix.T

def plot_with_landmarks(rcsb_id:str, mesh_grid_dimensions:np.ndarray, translation_vectors:np.ndarray, ):
    """
    @translation_vectors is a np.ndarray of shape (2,3) where 
        - the first row is the means of the coordinate set 
        - the second row is the deviations of the normalized coordinate set
        (to be used to reverse the normalization process or to travel to this coordinate frame)

    """
    with open( tunnel_atom_encoding_path(rcsb_id), "r", ) as infile:
        tunnel_atoms_data:list[dict] = json.load(infile)

    with open( ptc_data_path(rcsb_id), "r", ) as infile:
        ptc_data = json.load(infile)

    atom_coordinates_by_chain:dict[str,list] = { }
    for atom in tunnel_atoms_data:
        if atom['chain_nomenclature'][0] not in atom_coordinates_by_chain:
            atom_coordinates_by_chain[ atom['chain_nomenclature'][0] ] = []
        atom_coordinates_by_chain[ atom['chain_nomenclature'][0] ].extend( [ atom['coord'] ])

    ptc_midpoint = np.array(ptc_data["midpoint_coordinates"])
    print("GOT TRANSLATION VECTORS ", translation_vectors)
    # pprint(atom_coordinates_by_chain)
    # pprint(ptc_midpoint)
    #!----- convert to the normalized coord frame
    # ptc_midpoint = ptc_midpoint - translation_vectors[0] + translation_vectors[1]
    # for chain_name, coords in atom_coordinates_by_chain.items():
    #     atom_coordinates_by_chain[chain_name] = np.array(coords) - translation_vectors[0] + translation_vectors[1]
    # #!----- convert to the normalized coord frame

    mesh      = pv.read(poisson_recon_path(rcsb_id))
    plotter   = pv.Plotter()
    plotter.add_mesh(mesh, opacity=0.5)
    # plotter.show()
    colors = ['green','yellow', 'blue', 'magenta','cyan', 'pink', 'orange', 'purple', 'brown', 'grey']

    for i, ( chain_name, coords ) in enumerate(atom_coordinates_by_chain.items()):
        if chain_name == 'eL39':
            chain_color = 'blue'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=4, color=chain_color, render_points_as_spheres=True)
        if chain_name == 'uL4':
            chain_color =  'green'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=4, color=chain_color, render_points_as_spheres=True)
        if chain_name == 'uL22':
            chain_color =  'yellow'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=4, color=chain_color, render_points_as_spheres=True)
        else:
            continue

    plotter.add_points(move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array([ptc_midpoint])), point_size=6, color='red', render_points_as_spheres=True)
    plotter.add_text('Label Text', position='upper_left', font_size=18)
    plotter.show(auto_close=False)

def main():

    "the data arrives here as atom coordinates extracted from the biopython model "
    if not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)):
        bbox_atoms          = extract_bbox_atoms(RCSB_ID)
        bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms)
        print("Extracted tunnel atom encoding from PDB: {}.".format(RCSB_ID))
    else:
        with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile:
            bbox_atoms = json.load(infile)
        print("Opened tunnel atom encoding from file: {}.".format(tunnel_atom_encoding_path(RCSB_ID)))
        bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms)

    xyz_positive, xyz_negative, mesh_grid_dimensions, normalization_vectors = index_grid(bbox_atoms_expanded)
    np.save(normalization_vectors_path(RCSB_ID), normalization_vectors)
    db, clusters_container = interior_capture_DBSCAN(xyz_negative)

    surface_pts = surface_pts_via_convex_hull(clusters_container[DBSCAN_CLUSTER_ID])


    estimate_normals(surface_pts)
    plot_with_landmarks(RCSB_ID,mesh_grid_dimensions, normalization_vectors)


if __name__ == "__main__":
    main()
