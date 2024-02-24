from functools import partial
from pprint import pprint
from time import time
from matplotlib import pyplot as plt
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
