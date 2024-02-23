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
from mesh_generation.voxelize import expand_atomcenters_to_spheres_threadpool, get_sphere_indices_voxelized, normalize_atom_coordinates
from ribctl import  EXIT_TUNNEL_WORK

#? ---------- Params ------------

RCSB_ID = "6Z6K"

u_EPSILON     = 5.5
u_MIN_SAMPLES = 600
u_METRIC      = "euclidean"

DBSCAN_CLUSTER_ID = 3

#? ---------- Paths ------------

tunnel_atom_encoding_handle = "{}_tunnel_atoms_bbox.json".format(RCSB_ID)
tunnel_atom_encoding_path   = os.path.join(EXIT_TUNNEL_WORK,tunnel_atom_encoding_handle)
select_dbscan_cluster_path  = os.path.join(EXIT_TUNNEL_WORK, "dbscan_cluster_{}.npy".format(DBSCAN_CLUSTER_ID))

#? ------------------------------

def extract_bbox_atoms(rcsb_id:str)->list:

    centerline_expansion_atoms       = parse_struct_via_centerline(rcsb_id, open_tunnel_csv(rcsb_id))
    centerline_expansion_coordinates = np.array([a.get_coord() for a in centerline_expansion_atoms])
    bbox                             = bounding_box(centerline_expansion_coordinates)
    bbox_interior_atoms              = parse_struct_via_bbox(rcsb_id, bbox)

    return encode_atoms(rcsb_id, bbox_interior_atoms, write=True, writepath=tunnel_atom_encoding_path)

def expand_bbox_atoms_to_spheres(bbox_data:list|None ):

    if bbox_data is None:
        with open( tunnel_atom_encoding_path, "r", ) as infile:
            bbox_data = json.load(infile)

    __cords = np.array(list(map(lambda x: x["coord"], bbox_data)))
    __radii = np.array(list(map(lambda x: x["vdw_radius"], bbox_data)))

    sphere_sources = zip(__cords, __radii)

    SINK     = []
    expanded = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)


    return np.array(expanded)

def index_grid(expanded_sphere_voxels:np.ndarray, voxel_size:int=1):
    normalized_sphere_cords = normalize_atom_coordinates(expanded_sphere_voxels)
    sphere_cords_quantized  = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)

    max_values      = np.max(sphere_cords_quantized, axis=0)
    grid_dimensions = max_values + 1
    vox_grid        = np.zeros(grid_dimensions)

    vox_grid[
        sphere_cords_quantized[:, 0],
        sphere_cords_quantized[:, 1],
        sphere_cords_quantized[:, 2],
    ] = 1  

    xyz_v_negative = np.asarray( np.where(vox_grid != 1) )  
    xyz_v_positive = np.asarray( np.where(vox_grid == 1) )  # get back indexes of populated voxels

    return  xyz_v_positive,xyz_v_negative, vox_grid

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

    print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(xyz_v_negative.T), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) )
    db = DBSCAN(eps=u_EPSILON, min_samples=u_MIN_SAMPLES, metric=u_METRIC, n_jobs=5).fit( xyz_v_negative.T )
    labels = db.labels_
    dbscan_colors = [ cluster_colors[label * 2] if label != -1 else [0, 0, 0, 1] for label in labels ]

    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_    = list(labels).count(-1)
    print(" Estimated number of [ clusters, noise points ]: [ {}, {} ] ".format(n_clusters_, n_noise_ ) )

    CLUSTERS_CONTAINER = {}
    for point, label in zip(xyz_v_negative.T, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    for k, v in CLUSTERS_CONTAINER.items():
        print("Cluster {}: {} points".format(k, len(v)))

    return db, CLUSTERS_CONTAINER

if not os.path.exists(tunnel_atom_encoding_path):
    bbox_atoms          = extract_bbox_atoms(RCSB_ID)
    bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms)
else:
    with open( tunnel_atom_encoding_path, "r", ) as infile:
        bbox_atoms = json.load(infile)

    bbox_atoms_expanded = expand_bbox_atoms_to_spheres(bbox_atoms)

xyz_positive, xyz_negative, _ = index_grid(bbox_atoms_expanded)
db, clusters_container = interior_capture_DBSCAN(xyz_negative)


ptcloud_data_cluster = clusters_container[3]
ptcloud_data_positive = xyz_positive.T




def DBSCAN_CLUSTERS_visualize_all(dbscan_instance, dbscan_clusters):
    # create a dictionary of plt colors
    cluster_colors = dict(zip(range(-1, 40), plt.cm.terrain(np.linspace(0, 1, 40))))
    for k, v in cluster_colors.items():
        cluster_colors[k] = [*v[:3], 0.5]

    labels = dbscan_instance.labels_
    dbscan_colors = [ cluster_colors[label * 2] if label != -1 else [0, 0, 0, 1] for label in labels ]
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    n_noise_    = list(labels).count(-1)

    combined_clusters = ()
    point_cloud         = pv.PolyData()
    point_cloud["rgba"] = dbscan_clusters
    point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)


def DBSCAN_CLUSTERS_visualize_one(positive_space:np.ndarray, select_cluster:np.ndarray):
    rgbas_cluster  = [[15, 10, 221, 1] for datapoint in ptcloud_data_cluster]
    rgbas_positive = np.array([[205, 209, 228, 0.2] for _ in xyz_positive.T])

    combined       = np.concatenate([positive_space, select_cluster])
    rgbas_combined = np.concatenate([rgbas_cluster, rgbas_positive])

    point_cloud         = pv.PolyData(combined)
    point_cloud["rgba"] = rgbas_combined
    point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)

def estimate_normals():
    pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(ptcloud_data_cluster))
    pcd.estimate_normals( search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=40))


DBSCAN_CLUSTERS_visualize_all(db, clusters_container)



# def visualize_mesh():

    # # my_mesh = get_mesh(points3d,boundary_faces, opacity=1)
    # figb = go.Figure(data=[boundary_faces])
    # figb.update_layout(title_text="tit", title_x=0.5, width=800, height=800)
    # figb.show()