from pprint import pprint
import numpy as np

# from mesh_generation.visualization import visualize_pointcloud


# grid_dim = [10,20,40]
# grid=  np.zeros(grid_dim)


# some_points = np.array( [
#     [5,5,5],

#     [4,5,5],
#     [6,5,5],

#     [5,4,5],
#     [5,6,5],

#     [5,5,4],
#     [5,5,6],
# ] )


# grid[
#     some_points[:,0],
#     some_points[:,1],
#     some_points[:,2]
# ] = 1


# pprint(grid)


# positive = np.array(np.where(grid==1)).T
# negative = np.array(np.where(grid!=1)).T

# visualize_pointcloud(grid, "XXXX")
# visualize_pointcloud(positive, "XXXX")
# visualize_pointcloud(positive, "XXXX")








import pyvista as pv
def try_pyvista_pointcloud_axes():
    plotter               = pv.Plotter()
    # plotter.subplot(0,0)
    n_labels = 7
    # plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 8)

    # rgbas_cluster1       = [[15, 100, 21, 1] for datapoint in ptcloud1]
    # rgbas_positive1      = np.array([[205, 209, 228, 0.2] for _ in background_positive])

    combined1            = np.array([
        [100, 100, 100],
        [102,105,104],
        [90, 95,92]

    ])
    # rgbas_combined1      = np.concatenate([rgbas_cluster1, rgbas_positive1])

    point_cloud1         = pv.PolyData(combined1)
    # point_cloud1["rgba"] = rgbas_combined1

    plotter.add_mesh(point_cloud1)
    plotter.show()


try_pyvista_pointcloud_axes()