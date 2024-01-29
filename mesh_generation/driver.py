import open3d as o3d
import json
import os

import numpy as np
from ribctl.lib.tunnel import encode_atoms, open_tunnel_csv, parse_struct_via_centerline, get_sphere_indices_voxelized


RCSB_ID = '6Z6K'


tunnel_encoding = '/home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms.json'.format(RCSB_ID)
encoded         = None

if not os.path.isfile(tunnel_encoding):
    data    = open_tunnel_csv(RCSB_ID)
    atoms   = parse_struct_via_centerline(RCSB_ID,data, 15)
    encoded = encode_atoms(RCSB_ID,atoms)
    with open(tunnel_encoding, 'w') as f:
        json.dump(encoded,f)
        print("Wrote tunnel atoms to {}".format(f))

else:
    with open(tunnel_encoding, 'r') as f:
        encoded = json.load(f)

print(tunnel_encoding)
print("Processing {} tunnel atoms".format(len(encoded)))

spheres_sources = [(atom['coord'],atom['vdw_radius']*10) for atom in encoded]
expanded_coords = [[0,0,0]]
# expanded_coords = np.array([])


# for i,atom in enumerate(encoded):
#     if i % 100 == 0:
#         print("Processing atom {} of {}".format(i,len(encoded)))


# #!- --- -# 
# local_coordinates = get_sphere_indices_voxelized(coordinate, vdw_R*10)
# for _ in local_coordinates:
#     expanded_coords.append(_)
# #!- --- -#

import concurrent.futures
num_workers = 20

# Main function to coordinate the thread pool
def voxelize_to_spheres():
    # Function to be executed by each worker
    def sphere_task(args):
        coordinate, vdw_R = args
        result = get_sphere_indices_voxelized(coordinate, vdw_R)
        return result

    def update_expanded_coordinates(result):
        expanded_coords.extend(result)

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_args = {executor.submit(sphere_task, args): args for args in spheres_sources}
        for future in concurrent.futures.as_completed(future_to_args):
            result = future.result()
            print("Got {} results".format(len(result)))
            update_expanded_coordinates(result)




voxelize_to_spheres()



        # voxel_grid.add_voxel(o3d.geometry.Voxel(_), voxel_size=0.05)

pcd        = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(np.array(expanded_coords))
voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd, voxel_size=0.1)
o3d.visualization.draw_geometries([voxel_grid])



