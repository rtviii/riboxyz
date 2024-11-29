from functools import partial
from typing import Optional, Tuple
from matplotlib import pyplot as plt
import numpy as np
import concurrent.futures

from sklearn.neighbors import KDTree


def sphere_task(container_sink:list, atom_center_coordinate:np.ndarray, vdw_R=2):
    #TODO: Differentiate between atoms sizes
    result  = atompos_to_voxelized_sphere(atom_center_coordinate, 2)
    container_sink.extend(result)
    return result

def expand_atomcenters_to_spheres_threadpool(sink_container:list, atoms):
  with concurrent.futures.ThreadPoolExecutor(max_workers=12) as executor:
        futures = []
        for atom_pos in atoms:
            partial_task = partial(sphere_task, sink_container, atom_pos)
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

def atompos_to_voxelized_sphere(center: np.ndarray, radius: int):
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
    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint

    sphere_active_ix = []

    for ind in indices_list:
        x_ = ind[0]
        y_ = ind[1]
        z_ = ind[2]
        if (x_ - x0) ** 2 + (y_ - y0) ** 2 + (z_ - z0) ** 2 <= radius**2:
            sphere_active_ix.append([x_, y_, z_])

    return np.array(sphere_active_ix)

def index_grid(expanded_sphere_voxels: np.ndarray) :

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


import numpy as np
from sklearn.neighbors import KDTree

def make_cylinder_grid(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    resolution: float,
    existing_points: np.ndarray ,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Creates a cylindrical grid and tracks which points are "active" (occupied by existing points)
    """
    # Calculate number of samples for each dimension
    n_radial  = max(int(radius / resolution), 2)
    n_angular = max(int(2 * np.pi * radius / resolution), 8)
    n_height  = max(int(height / resolution), 2)
    
    # Create 1D coordinate arrays
    r = np.linspace(0, radius, n_radial)
    theta = np.linspace(0, 2*np.pi, n_angular)
    h = np.linspace(0, height, n_height)
    
    # Create 3D coordinate arrays
    R, Theta, H = np.meshgrid(r, theta, h, indexing='ij')
    
    # Convert to Cartesian coordinates
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)
    Z = H
    
    # Get the original shapes before flattening
    original_shape = X.shape
    
    # Stack and reshape coordinates into points array
    points = np.column_stack([
        X.ravel(), 
        Y.ravel(), 
        Z.ravel()
    ])
    
    # Transform points to align with cylinder axis
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    if axis_length == 0:
        raise ValueError("base_point and axis_point cannot be the same")
    
    axis = axis / axis_length
    
    # Create rotation matrix to align cylinder with given axis
    z_axis = np.array([0, 0, 1])
    if not np.allclose(axis, z_axis):
        rotation_axis = np.cross(z_axis, axis)
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        cos_theta = np.dot(z_axis, axis)
        sin_theta = np.sqrt(1 - cos_theta**2)
        
        K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                     [rotation_axis[2], 0, -rotation_axis[0]],
                     [-rotation_axis[1], rotation_axis[0], 0]])
        
        R = np.eye(3) + sin_theta * K + (1 - cos_theta) * (K @ K)
        points = points @ R.T
    
    # Translate points to base_point
    points = points + base_point
    
    # Initialize activity mask to match points shape
    activity_mask = np.zeros(len(points), dtype=bool)
    
    if existing_points is not None and len(existing_points) > 0:
        # Build KD-tree from existing points
        tree = KDTree(existing_points)
        # Query the tree for nearest neighbors
        distances, _ = tree.query(points, k=2)
        # Mark points as active if they're within resolution distance of existing points
        activity_mask = distances.ravel() <= resolution
    
    return points, activity_mask

def get_empty_space_points(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    resolution: float,
    existing_points: np.ndarray
) -> np.ndarray:
    """Get points representing empty space in the cylinder"""
    all_points, activity_mask = make_cylinder_grid(
        base_point, 
        axis_point, 
        radius, 
        height, 
        resolution, 
        existing_points
    )
    
    # Make sure mask shape matches points
    assert len(activity_mask) == len(all_points), \
        f"Mask shape {activity_mask.shape} doesn't match points shape {all_points.shape}"
    
    # Return points where there are no existing points nearby
    return all_points[~activity_mask]



# Example usage:
if __name__ == "__main__":
    # Example parameters
    base_point = np.array([0, 0, 0])
    axis_point = np.array([0, 0, 10])
    radius = 5
    height = 10
    resolution = 0.5
    
    # Generate some example existing points (random points in cylinder)
    n_points = 1000
    existing_points = np.random.rand(n_points, 3)
    existing_points[:, 2] *= height  # Scale z coordinates to height
    existing_points[:, :2] *= radius  # Scale x,y coordinates to radius
    
    # Get empty space points
    empty_points = get_empty_space_points(
        base_point, axis_point, radius, height,
        resolution, existing_points
    )
    
    print(f"Generated {len(empty_points)} empty space points")
    
    # Optional: Visualize results with matplotlib
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot existing points in blue
    ax.scatter(existing_points[:, 0], existing_points[:, 1], existing_points[:, 2], 
              c='blue', alpha=0.6, label='Existing Points')
    
    # Plot empty space points in red
    ax.scatter(empty_points[:, 0], empty_points[:, 1], empty_points[:, 2],
              c='red', alpha=0.2, label='Empty Space')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()