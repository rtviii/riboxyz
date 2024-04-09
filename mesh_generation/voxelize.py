from functools import partial
from matplotlib import pyplot as plt
import numpy as np
import concurrent.futures


def sphere_task(container_sink:list, atom_center_coordinate:np.ndarray, vdw_R=2):
    #TODO: Differentiate between atoms sizes
    result  = get_sphere_indices_voxelized(atom_center_coordinate, 2)
    container_sink.extend(result)
    return result

def expand_atomcenters_to_spheres_threadpool(sink_container:list, sphere_sources):
  with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        futures = []
        for (center_coordinate, vdw_R) in sphere_sources:
            partial_task = partial(sphere_task, sink_container, center_coordinate, vdw_R)
            future = executor.submit(partial_task)
            futures.append(future)
        concurrent.futures.wait(futures)

  return sink_container

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

    #!------ Generate indices of a voxel cube of side 2r  around the centerpoint
    x_range = slice(int(np.floor(x0 - radius)), int(np.ceil(x0 + radius)))
    y_range = slice(int(np.floor(y0 - radius)), int(np.ceil(y0 + radius)))
    z_range = slice(int(np.floor(z0 - radius)), int(np.ceil(z0 + radius)))

    indices = np.indices(
        (
            x_range.stop - x_range.start,
            y_range.stop - y_range.start,
            z_range.stop - z_range.start,
        )
    )

    indices      += np.array([x_range.start, y_range.start, z_range.start])[ :, np.newaxis, np.newaxis, np.newaxis ]
    indices       = indices.transpose(1, 2, 3, 0)
    indices_list  = list(map(tuple, indices.reshape(-1, 3)))
    print("Got indices list around ", center, " with radius ", radius, " and length: ", len(indices_list))
    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint

    sphere_active_ix = []

    for ind in indices_list:
        x_ = ind[0]
        y_ = ind[1]
        z_ = ind[2]
        if (x_ - x0) ** 2 + (y_ - y0) ** 2 + (z_ - z0) ** 2 <= radius**2:
            sphere_active_ix.append([x_, y_, z_])

    return np.array(sphere_active_ix)

def index_grid(expanded_sphere_voxels: np.ndarray, TRUNCATION_TUPLES:list[list|None]|None = None, normalize:bool = True) :

    def normalize_atom_coordinates(coordinates: np.ndarray)->tuple[ np.ndarray, np.ndarray ]:
        """@param coordinates: numpy array of shape (N,3)"""

        C      = coordinates
        mean_x = np.mean(C[:, 0])
        mean_y = np.mean(C[:, 1])
        mean_z = np.mean(C[:, 2])

        Cx = C[:, 0] - mean_x
        Cy = C[:, 1] - mean_y
        Cz = C[:, 2] - mean_z
        

        [dev_x, dev_y, dev_z] = [np.min(Cx), np.min(Cy), np.min(Cz)]

        #! shift to positive quadrant
        Cx = Cx + abs(dev_x)
        Cy = Cy + abs(dev_y)
        Cz = Cz + abs(dev_z)

        rescaled_coords = np.array(list(zip(Cx, Cy, Cz)))

        return rescaled_coords, np.array([[mean_x,mean_y,mean_z], [abs( dev_x ), abs( dev_y ), abs( dev_z )]])

    normalized_sphere_cords, mean_abs_vectors = normalize_atom_coordinates(expanded_sphere_voxels)
    voxel_size = 1

    sphere_cords_quantized = np.round(np.array(normalized_sphere_cords / voxel_size) ).astype(int)
    max_values             = np.max(sphere_cords_quantized, axis=0)
    grid_dimensions        = max_values + 1
    vox_grid               = np.zeros(grid_dimensions)

    print("Dimension of the voxel grid is ", vox_grid.shape)

    vox_grid[
        sphere_cords_quantized[:, 0],
        sphere_cords_quantized[:, 1],
        sphere_cords_quantized[:, 2]  ] = 1


    return ( vox_grid, grid_dimensions, mean_abs_vectors )


    if TRUNCATION_TUPLES is not None:
        if len(TRUNCATION_TUPLES) != 3:
            raise IndexError("You have to enter three truncation parameters for x, y, and z axes. (ex. ||20:50 to skip x and y axes )")

        if TRUNCATION_TUPLES[0] is not None:
            if TRUNCATION_TUPLES[0][0] is not None:
                vox_grid[:TRUNCATION_TUPLES[0][0],:,:]  = 1

            if TRUNCATION_TUPLES[0][1] is not None:
                 vox_grid[TRUNCATION_TUPLES[0][1]:,:,:] = 1

        if TRUNCATION_TUPLES[1] is not None:
            if TRUNCATION_TUPLES[1][0] is not None:
                vox_grid[:,:TRUNCATION_TUPLES[1][0],:]  = 1

            if TRUNCATION_TUPLES[1][1] is not None:
                 vox_grid[:,TRUNCATION_TUPLES[1][1]:,:] = 1

        if TRUNCATION_TUPLES[2] is not None:

            if TRUNCATION_TUPLES[2][0] is not None:
                vox_grid[:,:,:TRUNCATION_TUPLES[2][0]]  = 1

            if TRUNCATION_TUPLES[2][1] is not None:
                 vox_grid[:,:,TRUNCATION_TUPLES[2][1]:] = 1


    print("Voxel grid shape(post truncation if applied): ", vox_grid.shape)
    __xyz_v_negative_ix = np.asarray(np.where(vox_grid != 1))
    __xyz_v_positive_ix = np.asarray(np.where(vox_grid == 1))


    return (
        __xyz_v_positive_ix.T,
        __xyz_v_negative_ix.T,
        grid_dimensions,
        mean_abs_vectors )