import argparse
from pprint import pprint
import subprocess
import typing
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import json
import os
import numpy as np
import numpy as np
from mesh_generation.bbox_extraction import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_bbox,
    parse_struct_via_centerline,
)
from compas.geometry import bounding_box
from mesh_generation.paths import *
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
import warnings
warnings.filterwarnings("ignore")


def apply_poisson_reconstruction(surf_estimated_ptcloud_path: str, output_path: str, recon_depth:int=6, recon_pt_weight:int=3):
    import plyfile

    # command   = [ POISSON_RECON_BIN, "--in", surface_with_normals_path(rcsb_id), "--out", poisson_recon_path, "--depth", "6", "--pointWeight", "3", ]
    command = [
        POISSON_RECON_BIN,
        "--in",
        surf_estimated_ptcloud_path,
        "--out",
        output_path,
        "--depth",
        str(recon_depth),
        "--pointWeight",
        str(recon_pt_weight),
        # "--polygonMesh"
    ]
    process = subprocess.run(command, capture_output=True, text=True)

    if process.returncode == 0:
        print(">>PoissonRecon executed successfully.")
        print(">>Wrote {}".format(output_path))
        # Convert the plyfile to asciii
        data = plyfile.PlyData.read(output_path)
        data.text = True
        ascii_duplicate =output_path.split(".")[0] + "_ascii.ply"
        data.write(ascii_duplicate)
        print(">>Wrote {}".format(ascii_duplicate))
    else:
        print(">>Error:", process.stderr)

def ptcloud_convex_hull_points(pointcloud: np.ndarray, ALPHA:float, TOLERANCE:float) -> np.ndarray:
    assert pointcloud is not None

    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points

def estimate_normals(convex_hull_surface_pts: np.ndarray, output_path: str, kdtree_radius=None, kdtree_max_nn=None, correction_tangent_planes_n=None): 
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=kdtree_radius, max_nn=kdtree_max_nn) )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    o3d.io.write_point_cloud(output_path, pcd)
    return pcd
