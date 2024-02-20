import numpy as np
import pyvista as pv

# RCSB_ID = "6Z6K"
# # Load the PLY file
# file_path = "poisson_recon_{}.ply".format(RCSB_ID)  # Update with your file path

# mesh = pv.read(file_path)


# # Plot the mesh
# plotter = pv.Plotter()
# plotter.add_mesh(mesh)
# plotter.show()


#!----
# import open3d as o3d
# convex_hull = np.load("convex_hull_{}.npy".format(RCSB_ID))
# pcd        = o3d.geometry.PointCloud()
# pcd.points = o3d.utility.Vector3dVector(convex_hull)
# normals    = pcd.estimate_normals()

# o3d.visualization.draw_geometries([pcd])
# o3d.io.write_point_cloud("convex_hull_{}.ply".format(RCSB_ID), pcd)
# print("wrote", " convex_hull_{}.ply".format(RCSB_ID))



import open3d as o3d
import numpy as np

# Load the PLY file
ply_path = "convex_hull_6Z6K.ply"
point_cloud = o3d.io.read_point_cloud(ply_path)
o3d.visualization.draw_geometries([point_cloud], "Point Cloud Visualization")
point_cloud.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=15))
o3d.visualization.draw_geometries([point_cloud], 
                                   point_show_normal=True, 
                                   window_name="Point Cloud with Normals")
