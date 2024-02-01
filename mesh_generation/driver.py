from pprint import pprint
import open3d as o3d
import json
import os
import numpy as np
from ribctl.lib.tunnel import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_centerline,
    get_sphere_indices_voxelized,
)
import concurrent.futures

RCSB_ID = "6Z6K"


tunnel_encoding = (
    "/home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms.json".format(RCSB_ID)
)
encoded = None


def get_tunnel_atoms_file() -> list:
    if not os.path.isfile(tunnel_encoding):
        data = open_tunnel_csv(RCSB_ID)
        atoms = parse_struct_via_centerline(RCSB_ID, data, 15)
        encoded = encode_atoms(RCSB_ID, atoms)
        with open(tunnel_encoding, "w") as f:
            json.dump(encoded, f)
            print("Wrote tunnel atoms to {}".format(f))
            return encoded
    else:
        with open(tunnel_encoding, "r") as f:
            encoded = json.load(f)
            return encoded


# for i,atom in enumerate(encoded):
#     if i % 100 == 0:
#         print("Processing atom {} of {}".format(i,len(encoded)))


# #!- --- -#
# local_coordinates = get_sphere_indices_voxelized(coordinate, vdw_R*10)
# for _ in local_coordinates:
#     expanded_coords.append(_)
# #!- --- -#

num_workers = 14
encoded = get_tunnel_atoms_file()
# spheres_sources = [([ atom['coord'][0]*5, atom['coord'][1]*5, atom['coord'][2]*5 ], atom['vdw_radius']*5) for atom in encoded[:100]]
center_coordinates = np.array([atom["coord"] for atom in encoded])
spheres_sources = [(atom["coord"], atom["vdw_radius"]) for atom in encoded]
# expanded_coords = [[0, 0, 0]]
expanded_coords = []
print("Processing {} tunnel atoms".format(len(encoded)))


# def voxelize_to_spheres():
#     # Function to be executed by each worker
#     def sphere_task(args):
#         coordinate, vdw_R = args
#         result = get_sphere_indices_voxelized(coordinate, vdw_R)
#         return result

#     def update_expanded_coordinates(result):
#         expanded_coords.extend(result)

#     with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
#         future_to_args = {
#             executor.submit(sphere_task, args): args for args in spheres_sources
#         }

#         for future in concurrent.futures.as_completed(future_to_args):
#             result = future.result()
#             print(
#                 "Added {} coordinates to expanded radii. Total: {}".format(
#                     len(result), len(expanded_coords)
#                 )
#             )
#             update_expanded_coordinates(result)

# voxelize_to_spheres()

# np.array(expanded_coords)


pcd        = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(np.array(center_coordinates))
voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd, voxel_size=1)
pprint(dir(voxel_grid))

o3d.visualization.draw_geometries([pcd])
o3d.visualization.draw_geometries([voxel_grid])

import multiprocessing as mp
print(mp.cpu_count())

pprint(dir(o3d.geometry.VoxelGrid))
# import pyvista
# import numpy as np
# import pooch
# import pyvista as pv
# import PVGeo

# pv.global_theme.allow_empty_mesh = True
# # Create the 3D NumPy array of spatially referenced data
# # This is spatially referenced such that the grid is 20 by 5 by 10
# #   (nx by ny by nz)
# values = np.linspace(0, 10, 1000).reshape((20, 5, 10))
# print("Shape is",values.shape)
# values.shape

# # Create the spatial reference
# grid = pv.ImageData()

# # Set the grid dimensions: shape because we want to inject our values on the
# #   POINT data
# grid.dimensions = values.shape

# # Edit the spatial reference
# grid.origin = (0,0,0)  # The bottom left corner of the data set
# grid.spacing = (1, 1, 1)  # These are the cell sizes along each axis

# # Add the data values to the cell data
# grid.point_data["values"] = values.flatten(order="F")  # Flatten the array

# # Now plot the grid
# grid.plot(show_edges=True)


# pc = pv.PolyData(center_coordinates)

# vtkpoints = PVGeo.points_to_poly_data(center_coordinates)

# grid = PVGeo.filters.VoxelizePoints().apply(pc)

# p = pv.Plotter(notebook=0)
# p.add_mesh(grid, opacity=0.5, show_edges=True)
# p.add_mesh(pc, point_size=5, color="red")
# p.show_grid()

# p.show()
#! ----------------
# print(vtkpoints.points)
# vtkpoints.plot(clim=[0, 1], point_size=1)

# print(vtkpoints)
# voxelizer = PVGeo.filters.VoxelizePoints()

# grid = voxelizer.apply(vtkpoints)
# grid.plot()
# pc.glyph(geom=pv.Cube(), factor=1).plot()

# voxel_grid = pc.voxelize()

# p = pv.Plotter()
# print("Got plotter:" ,p)
# p.add_mesh(grid, opacity=0.5, show_edges=True)
# p.add_mesh(pc, point_size=5, color="red")
# p.show_grid()
# p.show()

# pdata                = pyvista.PolyData(np.array(center_coordinates))
# pdata["orig_sphere"] = np.arange(len(center_coordinates))

# # create many spheres from the point cloud
# sphere = pyvista.Sphere(radius=1.7, phi_resolution=10, theta_resolution=10)
# pc     = pdata.glyph(scale=False, geom=sphere, orient=False)
# pc.plot(cmap="Reds")

# voxel_grid.add_voxel(o3d.geometry.Voxel(_), voxel_size=0.05)
