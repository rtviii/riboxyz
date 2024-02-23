import argparse
from enum import Enum
from functools import partial
import json
import os
import pickle
from pprint import pprint
from time import time
import typing
import pyvista as pv
from matplotlib import pyplot as plt
import open3d as o3d
import sys
import numpy as np
import concurrent.futures






# Function to be executed by each worker
def sphere_task(container_sink:list, atom_center_coordinate:np.ndarray, vdw_R=2):
    result = get_sphere_indices_voxelized(atom_center_coordinate, 2)
    container_sink.extend(result)
    return result



def expand_atomcenters_to_spheres_threadpool(sink_container:list, sphere_sources):
  with concurrent.futures.ThreadPoolExecutor(max_workers=20) as executor:

        futures = []
        for (center_coordinate, vdw_R) in sphere_sources:
            partial_task = partial(sphere_task, sink_container, center_coordinate, vdw_R)
            future = executor.submit(partial_task)
            futures.append(future)

        concurrent.futures.wait(futures)

  return sink_container



DBSCAN_METRICS = [
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
    "cityblock",
    "cosine",
    "euclidean",
    "l1",
    "l2",
    "manhattan",
]

# parser = argparse.ArgumentParser(description="tunnel_extraction")

# parser.add_argument("--eps", type=float, help="Specify eps value")
# parser.add_argument("--min_samples", type=int, help="Specify min_samples value")
# parser.add_argument("--metric", choices=DBSCAN_METRICS, help="Choose a metric")

# parser.add_argument("--plot", action="store_true", help="Enable plotting")
# parser.add_argument("--cluster_only", type=int, help="Specify cluster only value")

# parser.add_argument("--rcsb_id", type=str, help="Specify RCSB ID", required=True)
# parser.add_argument("--input", type=str, help="Specify input file path")
# args = parser.parse_args()




def midpoints(x):
    sl = ()
    for _ in range(x.ndim):
        x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
        sl += np.index_exp[:]
    return x

def normalize_atom_coordinates(coordinates: np.ndarray):
    """@param coordinates: numpy array of shape (N,3)"""

    C = coordinates
    Cx = C[:, 0] - np.mean(C[:, 0])
    Cy = C[:, 1] - np.mean(C[:, 1])
    Cz = C[:, 2] - np.mean(C[:, 2])

    [dev_x, dev_y, dev_z] = [np.min(Cx), np.min(Cy), np.min(Cz)]

    #! shift to positive quadrant
    Cx = Cx + abs(dev_x)
    Cy = Cy + abs(dev_y)
    Cz = Cz + abs(dev_z)

    rescaled_coords = np.array(list(zip(Cx, Cy, Cz)))

    return rescaled_coords

def visualize_source_coordinates(
    nulled_grid: np.ndarray,
    coordinates: np.ndarray,
):
    for coordinate in coordinates:
        # coordinates of the side of the given voxel
        vox_x, vox_y, vox_z = (
            int(np.floor(coordinate[0])),
            int(np.floor(coordinate[1])),
            int(np.floor(coordinate[2])),
        )
        nulled_grid[vox_x, vox_y, vox_z] = True
    return nulled_grid

def plt_plot(x_ix, y_ix, z_ix, filled_grid):
    # facecolors =  np.zeros(filled.shape + (3,))
    ax = plt.figure().add_subplot(projection="3d")
    ax.voxels(x_ix, y_ix, z_ix, filled_grid, facecolors=[0, 1, 1, 0.3], linewidth=0.5)
    ax.set(xlabel="r", ylabel="g", zlabel="b")
    ax.set_aspect("equal")

    plt.show()
    exit()

def get_sphere_indices_voxelized(center: np.ndarray, radius: int):
    """Make sure radius reflects the size of the underlying voxel grid"""
    x0, y0, z0 = center

    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint
    x_range = slice(int(np.floor(x0) - (radius)), int(np.ceil(x0) + (radius)))
    y_range = slice(int(np.floor(y0) - (radius)), int(np.ceil(y0) + (radius)))
    z_range = slice(int(np.floor(z0) - (radius)), int(np.ceil(z0) + (radius)))

    indices = np.indices(
        (
            x_range.stop - x_range.start,
            y_range.stop - y_range.start,
            z_range.stop - z_range.start,
        )
    )

    indices += np.array([x_range.start, y_range.start, z_range.start])[
        :, np.newaxis, np.newaxis, np.newaxis
    ]
    indices = indices.transpose(1, 2, 3, 0)
    indices_list = list(map(tuple, indices.reshape(-1, 3)))
    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint

    sphere_active_ix = []

    for ind in indices_list:
        x_ = ind[0]
        y_ = ind[1]
        z_ = ind[2]
        if (x_ - x0) ** 2 + (y_ - y0) ** 2 + (z_ - z0) ** 2 <= radius**2:
            sphere_active_ix.append([x_, y_, z_])

    return np.array(sphere_active_ix)

# cords_SPHERES = []

# def voxelize_to_spheres(sphere_sources):

#     num_workers = 20

#     # Function to be executed by each worker
#     def sphere_task(args):
#         coordinate, vdw_R = args
#         result = get_sphere_indices_voxelized(coordinate, 2)
#         return result

#     def update_expanded_coordinates(result):
#         cords_SPHERES.extend(result)

#     with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
#         future_to_args = {
#             executor.submit(sphere_task, args): args for args in sphere_sources
#         }

#         for future in concurrent.futures.as_completed(future_to_args):
#             result = future.result()
#             update_expanded_coordinates(result)


# RCSB_ID = ""
# with open( "/home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms_bbox.json".format( RCSB_ID ), "r", ) as infile:
#     ptcloud_data_cluster = json.load(infile)

# __cords = np.array(list(map(lambda x: x["coord"], ptcloud_data_cluster)))
# __radii = np.array(list(map(lambda x: x["vdw_radius"], ptcloud_data_cluster)))

# sphere_sources = zip(__cords, __radii)

# sink = []
# expanded = expand_atomcenters_to_spheres_threadpool(sink, sphere_sources)


# normalized_sphere_cords = normalize_atom_coordinates(cords_SPHERES)
# voxel_size              = 1
# sphere_cords_quantized  = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)

# max_values      = np.max(sphere_cords_quantized, axis=0)
# grid_dimensions = max_values + 1
# vox_grid        = np.zeros(grid_dimensions)

# vox_grid[
#     sphere_cords_quantized[:, 0],
#     sphere_cords_quantized[:, 1],
#     sphere_cords_quantized[:, 2],
# ] = 1  

# xyz_v_negative = np.asarray( np.where(vox_grid != 1) )  
# xyz_v_positive = np.asarray( np.where(vox_grid == 1) )  # get back indexes of populated voxels


#! --------------------------------- DBSCAN


# From scipy.spatial.distance:
metrics = [
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis",
    "minkowski",
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean",
    "yule",
    "cityblock",
    "cosine",
    "euclidean",
    "l1",
    "l2",
    "manhattan",
]

cluster_colors = dict(zip(range(-1, 40), plt.cm.terrain(np.linspace(0, 1, 40))))
for k, v in cluster_colors.items():
    cluster_colors[k] = [*v[:3], 0.5]

import numpy as np
from sklearn.cluster import DBSCAN


# attempts: list[typing.Tuple[float, int, str]] = [
#     (2.5, 30, "sqeuclidean"),
#     (
#         2,
#         20,
#         "euclidean",
#     ),  # In 0.407s. Estimated number of [ clusters, noise points ]: [ 22, 6389 ]
#     (2, 40, "euclidean"),
#     (
#         3,
#         40,
#         "euclidean",
#     ),  # In 0.602s. Estimated number of [ clusters, noise points ]: [ 3, 1040 ]
#     (2, 60, "euclidean"),
#     (
#         4,
#         100,
#         "euclidean",
#     ),  # In 0.84s. Estimated number of [ clusters, noise points ]: [ 3, 2610 ]
#     (
#         5,
#         200,
#         "euclidean",
#     ),  # In 0.84s. Estimated number of [ clusters, noise points ]: [ 3, 2610 ]
#     (
#         6,
#         400,
#         "euclidean",
#     ),  # In 1.583s. Estimated number of [ clusters, noise points ]: [ 2, 8442 ]
#     (
#         5,
#         500,
#         "euclidean",
#     ),  # <<<< Really good result In 0.98s. Estimated number of [ clusters, noise points ]: [ 10, 123078 ]
#     (
#         6,
#         600,
#         "euclidean",
#     ),  # <<<< Best so far In 1.46s. Estimated number of [ clusters, noise points ]: [ 9, 53038 ]
#     (5.5, 600, "euclidean"),  # <<<< A little finer res
#     (5.5, 650, "euclidean"),  #
#     (5.5, 650, "canberra"),  # failed
#     (5.5, 650, "sqeuclidean"),  # 0 clusters
#     (4, 100, "sqeuclidean"),
#     (3, 40, "canberra"),
#     ( 6, 600, "euclidean", ),  # <<<< Best so far In 1.46s. Estimated number of [ clusters, noise points ]: [ 9, 53038 ]
#     (5.5, 600, "euclidean"),  # <<<< A little finer res
# ]

# u_EPSILON     = attempts[-1][0]
# u_MIN_SAMPLES = attempts[-1][1]
# u_METRIC      = attempts[-1][2]


# print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(xyz_v_negative.T), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) )
# db = DBSCAN(eps=u_EPSILON, min_samples=u_MIN_SAMPLES, metric=u_METRIC, n_jobs=5).fit(
#     xyz_v_negative.T
# )
# labels = db.labels_
# dbscan_colors = [ cluster_colors[label * 2] if label != -1 else [0, 0, 0, 1] for label in labels ]
# # Number of clusters in labels, ignoring noise if present.
# n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
# n_noise_    = list(labels).count(-1)
# print( "In {}s. Estimated number of [ clusters, noise points ]: [ {}, {} ] ".format( np.round(t2 - t1, 3), n_clusters_, n_noise_ ) )

# CLUSTER_CONTAINER = {}
# for point, label in zip(xyz_v_negative.T, labels):
#     if label not in CLUSTER_CONTAINER:
#         CLUSTER_CONTAINER[label] = []
#     CLUSTER_CONTAINER[label].append(point)

# for k, v in CLUSTER_CONTAINER.items():
#     print("Cluster {}: {} points".format(k, len(v)))

# exit()
# # Now to partition the points into clusters and noise


# if args.plot == True:
#     if args.cluster_only != None:
#         ptcloud_data_cluster = CLUSTER_CONTAINER[args.cluster_only]
#         ptcloud_data_positive = np.array(xyz_v_positive.T)


#         rgbas_cluster = [[15, 10, 221, 1] for datapoint in ptcloud_data_cluster]
#         rgbas_positive = np.array([[205, 209, 228, 0.2] for _ in xyz_v_positive.T])
#         # rgba_combined         = np.concatenate([rgba_p, rgba_n])
#         # point_cloud['rgba']  = rgba_combined
#         combined = np.concatenate([ptcloud_data_cluster, ptcloud_data_positive])
#         rgbas_combined = np.concatenate([rgbas_cluster, rgbas_positive])

#         point_cloud         = pv.PolyData(combined)
#         point_cloud["rgba"] = rgbas_combined
#         # point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)



#         cluster_path = "/home/rtviii/dev/riboxyz/mesh_generation/{}_cluster.npy".format(RCSB_ID)
#         np.save(cluster_path, ptcloud_data_cluster)
#         # !-----------

#         # pcd        = o3d.geometry.PointCloud()
#         # pcd.points = o3d.utility.Vector3dVector(ptcloud_data_cluster)
#         # normals    = pcd.estimate_normals()

#         o3d.visualization.draw_geometries([normals])
#         o3d.io.write_point_cloud("{}_cluster_cloud.ply".format(RCSB_ID), pcd)
#         print("wrote")
#         # cloud.plot(point_size=1)

#         # exit()
#         # #? Alpha of about >  2.5 starts to hide the detail
        # cloud = pv.PolyData(ptcloud_data_cluster)
        # grid              = cloud.delaunay_3d(alpha=3, tol=1.5,offset=2,progress_bar=True)
        # convex_hull = grid.extract_surface().cast_to_pointset()
        # print(convex_hull.points)
        # print(convex_hull.points.shape)
        # from ribctl import ASSETS_PATH
        # np.save("convex_hull_{}.npy".format(RCSB_ID), convex_hull.points)
        # print("Saved generated convex hull to ")
        # convex_hull.plot(show_edges=True)

        
        # grid.plot(show_edges=True)
        # edges = grid.extract_all_edges()
#         exit()

#         print(edges)
#         print(edges.points.shape)
#         print(edges.points)
#         exit()

#         surface = surf.extract_surface()
#         print(surface)
#         mask = surface.select_enclosed_points(ptcloud_data_cluster)
#         masked_points = cloud.points[mask]
#         print(cloud.points.shape)
#         print(masked_points.shape)

#         # from plotly import graph_objects as go
#         # boundary_mesh = surf.extract_geometry()
#         # boundary_faces = boundary_mesh.faces.reshape((-1,4))[:, 1:]  
#         # print("boundary faces shape: ",boundary_faces.shape)
#         # np.save("boundary_faces.npy", boundary_faces)
#         # #? Alpha of about >  2.5 starts to hide the detail


#         boundary_faces = np.load("boundary_faces.npy")
#         plotter = pv.Plotter()
#         plotter.add_points(boundary_faces, color='red', point_size=2)
#         plotter.show()



# #         # print(boundary_faces)
# #         # # my_mesh = get_mesh(points3d,boundary_faces, opacity=1)
# #         # figb = go.Figure(data=[boundary_faces])
# #         # figb.update_layout(title_text="tit", title_x=0.5, width=800, height=800)
# #         # figb.show()

# #         exit()

# #         pprint(surf)
# #         surf.save("{}_tunnel_delaunay.vtk".format(RCSB_ID))
# #         print(surf)
# #         surf.plot(show_edges=True)

# #         pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(ptcloud_data_cluster))
# #         pcd.estimate_normals( search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=40))

# #         # # o3d.visualization.draw_geometries([pcd],
# #         # #                                   zoom=0.3412,
# #         # #                                   front=[0.4257, -0.2125, -0.8795],
# #         # #                                   lookat=[2.6172, 2.0475, 1.532],
# #         # #                                   up=[-0.0694, -0.9768, 0.2024],
# #         # #                                   point_show_normal=True)
# #         # # exit()
# #         # def poisson_mesh(pc, depth=8, width=0, scale=1.1, linear_fit=False):
# #         #     """`pc` is a `pyvista.PolyData` point cloud. The default arguments are abitrary"""
# #         #     cloud = o3d.geometry.PointCloud()
# #         #     cloud.points = o3d.utility.Vector3dVector(pc.points)
# #         #     cloud.normals = o3d.utility.Vector3dVector(pc.normals)
# #         #     trimesh, _ = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(cloud, depth=depth, width=width, scale=scale, linear_fit=linear_fit)
# #         #     v = np.asarray(trimesh.vertices)
# #         #     f = np.array(trimesh.triangles)
# #         #     f = np.c_[np.full(len(f), 3), f]
# #         #     mesh = pv.PolyData(v, f)
# #         #     return mesh.clean()

# #         # x = poisson_mesh(pcd, depth=5, width=1, scale=1.1, linear_fit=False)
# #         # x.plot()
# #         # # print(x)

# #     else:
# #         # pv.set_plot_theme('dark')
# #         # ptcloud_data_combined = np.concatenate([xyz_v_positive.T, xyz_v_negative.T])
# #         ptcloud_data_negative = np.concatenate(xyz_v_negative.T)
# #         point_cloud = pv.PolyData(ptcloud_data_negative)
# #         # rgba_n                = np.array([[250,0,0,0.1] for _ in xyz_v_negative.T] )
# #         # rgba_p                = np.array([[10,200,200, 1] for _ in xyz_v_positive.T] )
# #         # rgba_combined         = np.concatenate([rgba_p, rgba_n])
# #         # point_cloud['rgba']  = rgba_combined

# #         point_cloud["rgba"] = dbscan_colors
# #         point_cloud.plot(scalars="rgba", rgb=True, notebook=False, show_bounds=True)


