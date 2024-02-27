import argparse
from pprint import pprint
import subprocess
import sys
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
from ribctl import  EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA



# DBSCAN_CLUSTER_ID = 3

DBSCAN_METRICS = [ "braycurtis", "canberra", "chebyshev", "correlation", "dice", 
 "hamming", "jaccard", "kulsinski", "mahalanobis", "minkowski",
  "rogerstanimoto", "russellrao", "seuclidean", "sokalmichener", 
  "sokalsneath", "sqeuclidean", "yule","cityblock", "cosine", "euclidean", "l1", "l2", "manhattan" ]
#? ---------- Paths ------------
tunnel_atom_encoding_path    = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_tunnel_atoms_bbox.json".format(rcsb_id.upper()))
normalization_vectors_path   = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_normalization_vectors.npy".format(rcsb_id.upper()))
selected_dbscan_cluster_path = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_dbscan_cluster.npy".format(rcsb_id.upper()))
convex_hull_cluster_path     = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_convex_hull.npy".format(rcsb_id.upper()))
surface_with_normals_path    = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_normal_estimated_surf.ply".format(rcsb_id.upper()))
poisson_recon_path           = lambda rcsb_id : os.path.join(EXIT_TUNNEL_WORK,    rcsb_id.upper(),"{}_poisson_recon.ply".format(rcsb_id.upper()))
ptc_data_path                = lambda rcsb_id : os.path.join(RIBETL_DATA,rcsb_id, rcsb_id.upper(),"{}_PTC_COORDINATES.json".format(rcsb_id.upper()))
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

    return encode_atoms(rcsb_id, bbox_interior_atoms, write=True, writepath=tunnel_atom_encoding_path(rcsb_id))

def expand_bbox_atoms_to_spheres(bbox_data:list|None, rcsb_id:str ):

    if bbox_data is None:
        with open( tunnel_atom_encoding_path(rcsb_id), "r", ) as infile:
            bbox_data = json.load(infile)

    __cords = np.array(list(map(lambda x: x["coord"], bbox_data)))
    __radii = np.array(list(map(lambda x: x["vdw_radius"], bbox_data)))

    sphere_sources = zip(__cords, __radii)
    SINK     = []
    expanded = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)

    return np.array(expanded)

def index_grid(expanded_sphere_voxels:np.ndarray):


    normalized_sphere_cords, mean_abs_vectors = normalize_atom_coordinates(expanded_sphere_voxels)

    voxel_size             = 1
    sphere_cords_quantized = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)

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

def interior_capture_DBSCAN(xyz_v_negative:np.ndarray, eps:float=5.5, min_samples:int=600, metric:str="euclidean"):

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
    db                 = DBSCAN(eps=eps, min_samples=min_samples, metric=metric, n_jobs=5).fit( xyz_v_negative )
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

def surface_pts_via_convex_hull(rcsb_id:str,selected_cluster:np.ndarray|None=None):

    if not os.path.exists(convex_hull_cluster_path(rcsb_id)):
        assert(selected_cluster is not None)
        cloud       = pv.PolyData(selected_cluster)
        grid        = cloud.delaunay_3d(alpha=2, tol=1.5,offset=2,progress_bar=True)
        convex_hull = grid.extract_surface().cast_to_pointset()
        np.save(convex_hull_cluster_path(rcsb_id), convex_hull.points)
        print("Saved generated convex hull to {}".format(convex_hull_cluster_path(rcsb_id)))
        return convex_hull.points
    else:
        pts = np.load(convex_hull_cluster_path(rcsb_id))
        print("Loaded generated convex hull from {}".format(convex_hull_cluster_path(rcsb_id)))
        return pts
    # convex_hull.plot(show_edges=True)

def estimate_normals(rcsb_id,convex_hull_surface_pts:np.ndarray|None=None):
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=15))
    pcd.orient_normals_consistent_tangent_plane(k=10)
    o3d.visualization.draw_geometries([pcd], point_show_normal=True, window_name="Point Cloud with Normals")
    o3d.io.write_point_cloud(surface_with_normals_path(rcsb_id), pcd)
    print("Wrote surface with normals {}".format(surface_with_normals_path(rcsb_id)))

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

    mesh      = pv.read(poisson_recon_path(rcsb_id))
    plotter   = pv.Plotter()
    plotter.add_mesh(mesh, opacity=0.5)

    CHAIN_PT_SIZE = 8
    PTC_PT_SIZE   = 20

    for i, ( chain_name, coords ) in enumerate(atom_coordinates_by_chain.items()):
        if chain_name == 'eL39':
            chain_color = 'blue'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=CHAIN_PT_SIZE, color=chain_color, render_points_as_spheres=True)
        if chain_name == 'uL4':
            chain_color =  'green'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=CHAIN_PT_SIZE, color=chain_color, render_points_as_spheres=True)
        if chain_name == 'uL22':
            chain_color =  'yellow'
            plotter.add_points( move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array(coords)),  point_size=CHAIN_PT_SIZE, color=chain_color, render_points_as_spheres=True)
        else:
            continue

    plotter.add_points(move_cords_to_normalized_cord_frame(mesh_grid_dimensions, translation_vectors, np.array([ptc_midpoint])), point_size=PTC_PT_SIZE, color='red', render_points_as_spheres=True)
    labels_colors = [('uL4', 'green'),('uL22','yellow'),('eL39','blue'), ('PTC','red')]

    for ( i, ( label, color ) ) in enumerate(labels_colors):
        offset = i * 30  # Adjust the offset as needed
        position = (20, 200 - offset, 0)
        plotter.add_text(label,
                     position             = position,
                     font_size            = 20,
                     color                = color,
                     shadow               = True)
    plotter.show(auto_close=False)

def apply_poisson_reconstruction(rcsb_id:str):
    command = [
        POISSON_RECON_BIN,
        "--in",
        surface_with_normals_path(rcsb_id),
        "--out",
        poisson_recon_path(rcsb_id),
        "--depth",
        "6",
        "--pointWeight",
        "3"
    ]

    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        print("PoissonRecon executed successfully.")
        print("Saved to {}".format(poisson_recon_path(rcsb_id)))
    else:
        print("Error:", process.stderr)


def main():

    #? ---------- Params ------------
    parser = argparse.ArgumentParser()
    # Add command-line arguments
    parser.add_argument("--rcsb_id"    , type   =str           , help="Specify the value for eps (float)"        , required=True)
    parser.add_argument("--eps"        , type   =float         , help="Specify the value for eps (float)"        )
    parser.add_argument("--min_samples", type   =int           , help="Specify the value for min_samples (int)"  )
    parser.add_argument("--metric"     , choices=DBSCAN_METRICS, help="Choose a metric from the provided options")

    # Parse the command-line arguments
    args = parser.parse_args()


    _u_EPSILON     = 5.5; _u_MIN_SAMPLES = 600; _u_METRIC      = "euclidean"

    RCSB_ID       = args.rcsb_id.upper()
    u_EPSILON     = args.eps if args.eps is not None else _u_EPSILON
    u_MIN_SAMPLES = args.min_samples if args.min_samples is not None else _u_MIN_SAMPLES
    u_METRIC      = args.metric if args.metric is not None else _u_METRIC

    struct_tunnel_dir = os.path.join(EXIT_TUNNEL_WORK, RCSB_ID)
    if not os.path.exists(struct_tunnel_dir):
         os.mkdir(struct_tunnel_dir)

    "the data arrives here as atom coordinates extracted from the biopython model "
    if not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)):
        bbox_atoms          = extract_bbox_atoms(RCSB_ID)
        bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms, RCSB_ID)
        print("Extracted tunnel atom encoding from PDB: {}.".format(RCSB_ID))
    else:
        with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile:
            bbox_atoms = json.load(infile)
        print("Opened tunnel atom encoding from file: {}.".format(tunnel_atom_encoding_path(RCSB_ID)))
        bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms, RCSB_ID)

    xyz_positive, xyz_negative, mesh_grid_dimensions, normalization_vectors = index_grid(bbox_atoms_expanded)
    np.save(normalization_vectors_path(RCSB_ID), normalization_vectors)
    db, clusters_container = interior_capture_DBSCAN(xyz_negative, u_EPSILON, u_MIN_SAMPLES, u_METRIC)



    DBSCAN_CLUSTERS_visualize_all(clusters_container)


    while True:
        user_input = input("Enter 'Q' to quit or a number between 1 and 20: ")
        if user_input.upper() == 'Q':
            break  # Exit the loop if the user enters 'Q'
        elif user_input.isdigit() and 0 <= int(user_input) <= 20:
            DBSCAN_CLUSTER_ID = int(user_input)
            surface_pts       = surface_pts_via_convex_hull(RCSB_ID,clusters_container[DBSCAN_CLUSTER_ID])

            estimate_normals(RCSB_ID,surface_pts)
            apply_poisson_reconstruction(RCSB_ID)
            plot_with_landmarks(RCSB_ID,mesh_grid_dimensions, normalization_vectors)
        else:
            print("Invalid input. Please enter 'Q' or a number between 1 and 20.")
    




if __name__ == "__main__":
    main()
