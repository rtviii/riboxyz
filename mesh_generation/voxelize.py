import json
import open3d as o3d
import sys
import numpy as np

RCSB_ID = sys.argv[1].upper()

with open( "/home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms_bbox.json".format(RCSB_ID), "r" ) as infile:
    data = json.load(infile)


num_workers = 20
__cords          = np.array(list(map(lambda x: x["coord"], data)))



pcd = o3d.geometry.PointCloud()
pcd.points = o3d.utility.Vector3dVector(__cords)
voxel_grid = o3d.geometry.VoxelGrid.create_from_point_cloud(pcd, voxel_size=1)
o3d.visualization.draw_geometries([pcd])