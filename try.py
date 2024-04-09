from pprint import pprint
import numpy as np

from mesh_generation.visualization import visualize_pointcloud


grid_dim = [10,20,40]
grid=  np.zeros(grid_dim)


some_points = np.array( [
    [5,5,5],

    [4,5,5],
    [6,5,5],

    [5,4,5],
    [5,6,5],

    [5,5,4],
    [5,5,6],
] )


grid[
    some_points[:,0],
    some_points[:,1],
    some_points[:,2]
] = 1


pprint(grid)


positive = np.array(np.where(grid==1)).T
negative = np.array(np.where(grid!=1)).T

# visualize_pointcloud(grid, "XXXX")
# visualize_pointcloud(positive, "XXXX")
# visualize_pointcloud(positive, "XXXX")






