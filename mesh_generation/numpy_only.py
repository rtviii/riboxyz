import json
import concurrent.futures
import numpy as np


def midpoints(x):
    sl = ()
    for _ in range(x.ndim):
        x = (x[sl + np.index_exp[:-1]] + x[sl + np.index_exp[1:]]) / 2.0
        sl += np.index_exp[:]
    return x


def get_sphere_indices_voxelized(center: np.ndarray, radius: int):
    """Make sure radius reflects the size of the underlying voxel grid"""
    x0, y0, z0 = center

    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint
    x_range = slice(int(np.floor(x0) - (radius + 1)), int(np.ceil(x0) + (radius + 1)))
    y_range = slice(int(np.floor(y0) - (radius + 1)), int(np.ceil(y0) + (radius + 1)))
    z_range = slice(int(np.floor(z0) - (radius + 1)), int(np.ceil(z0) + (radius + 1)))

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


def normalize_atom_coordinates(coordinates: np.ndarray):
    """@param coordinates: numpy array of shape (N,3) for atoms lining the tunnel in radius R (usually 15Angstrom) from MOLE centerline"""
    #! normalize to origin
    C = coordinates
    Cx = C[:, 0] - np.mean(C[:, 0])
    Cy = C[:, 1] - np.mean(C[:, 1])
    Cz = C[:, 2] - np.mean(C[:, 2])

    #! negative deviation from zero"
    dev = np.min([np.min(Cx), np.min(Cy), np.min(Cz)])

    #! shift to positive quadrant
    Cx = Cx + abs(dev)
    Cy = Cy + abs(dev)
    Cz = Cz + abs(dev)
    rescaled_coords = np.array(list(zip(Cx, Cy, Cz)))

    # ! Create a 3D grid 10 voxels larger than the max amplitude of the point cloud in any direction
    # amplitude_X = np.max(Cx) - np.min(Cx)
    # amplitude_Y = np.max(Cy) - np.min(Cy)
    # amplitude_Z = np.max(Cz) - np.min(Cz)

    # biggest_dimension = int( np.ceil(np.max([amplitude_Z, amplitude_Y, amplitude_X])) + 1 )
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


with open( "/home/rtviii/dev/riboxyz/mesh_generation/6Z6K_tunnel_atoms.json", "r" ) as infile:
    data = json.load(infile)


num_workers = 20

__cords          = np.array(list(map(lambda x: x["coord"], data)))
__radii          = np.array(list(map(lambda x: x["vdw_radius"], data)))
cords_NORMALIZED = normalize_atom_coordinates(__cords)
cords_SPHERES    = []
sphere_sources   = zip(cords_NORMALIZED, __radii)


def voxelize_to_spheres():
    # Function to be executed by each worker
    def sphere_task(args):
        coordinate, vdw_R = args
        result = get_sphere_indices_voxelized(coordinate, vdw_R)
        return result

    def update_expanded_coordinates(result):
        cords_SPHERES.extend(result)

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_args = {
            executor.submit(sphere_task, args): args for args in sphere_sources
        }

        for future in concurrent.futures.as_completed(future_to_args):
            result = future.result()
            print(
                "Added {} coordinates to expanded radii. Total: {}".format(
                    len(result), len(cords_SPHERES)
                )
            )
            update_expanded_coordinates(result)


voxelize_to_spheres()

cords_SPHERES = np.array(cords_SPHERES)

VOXEL_SIZE    = 1

coordinates_x = cords_SPHERES[:, 0]
max_dim_x     = np.max(coordinates_x)
coordinates_y = cords_SPHERES[:, 1]
max_dim_y     = np.max(coordinates_y)
coordinates_z = cords_SPHERES[:, 2]
max_dim_z     = np.max(coordinates_z)

xyz_q = np.round(np.array(cords_SPHERES / VOXEL_SIZE)).astype(
    int
)  # quantized point values, here you will loose precision
vox_grid = np.zeros(
    (
        int(max_dim_x / VOXEL_SIZE) + 1,
        int(max_dim_y / VOXEL_SIZE) + 1,
        int(max_dim_z / VOXEL_SIZE) + 1,
    )
)

vox_grid[xyz_q[:, 0], xyz_q[:, 1], xyz_q[:, 2]] = 1

xyz_v = np.asarray(np.where(vox_grid != 1))
final_source = xyz_v.T

import pyvista as pv

ptcloud = pv.PolyData(xyz_v.T)
plotter = pv.Plotter()

plotter.add_points(ptcloud, opacity=0.1)
plotter.show()

# exit()


# import open3d as o3d
# pcd = o3d.geometry.PointCloud()
# pcd.points = o3d.utility.Vector3dVector(xyz_v.T)
# pcd.colors = o3d.utility.Vector3dVector()
# voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd, voxel_size=1)
# o3d.visualization.draw_geometries([pcd])
# o3d.visualization.draw_geometries([voxel_grid])


# -------------;


# voxel_size  = 0.5

# rescaled_coordinates, dim = normalize_atom_coordinates(coordinates)


# x,y,z = np.indices((dim, dim, dim))
# xc = midpoints(x)
# yc = midpoints(y)
# zc = midpoints(z)
# filled = xc + yc + zc  < -1 # nulling out the entire grid

# for coordinate in cords_NORMALIZED:
#     # coordinates of the side of the given voxel
#     vox_x, vox_y, vox_z = (
#         int(np.floor(coordinate[0])),
#         int(np.floor(coordinate[1])),
#         int(np.floor(coordinate[2])),
#     )
#     nulled_grid[vox_x, vox_y, vox_z] = True
# filled = visualize_source_coordinates(filled,rescaled_coordinates)

# print(np.shape(filled))
