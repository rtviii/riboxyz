import argparse
from pprint import pprint
import subprocess
import typing
from matplotlib import pyplot as plt
import open3d as o3d
from Bio.PDB import MMCIFParser
from Bio.PDB import NeighborSearch
import pyvista as pv
import json
import os
import numpy as np
import numpy as np
from sklearn.cluster import DBSCAN
from __archive.scripts.pymol_visualtion import extract_chains
from mesh_generation.bbox_extraction import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_bbox,
    parse_struct_via_centerline,
)
from compas.geometry import bounding_box
from mesh_generation.visualization import (
    DBSCAN_CLUSTERS_visualize_largest,
    custom_cluster_recon_path,
    plot_multiple_by_kingdom,
    plot_multiple_surfaces,
    plot_with_landmarks,
    DBSCAN_CLUSTERS_particular_eps_minnbrs,
)
from mesh_generation.paths import *
from mesh_generation.voxelize import (
    expand_atomcenters_to_spheres_threadpool,
    normalize_atom_coordinates,
)
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_chains_by_auth_asym_id


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
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        print("PoissonRecon executed successfully.")
        print("Wrote {}".format(output_path))
        # Convert the plyfile to asciii
        data = plyfile.PlyData.read(surf_estimated_ptcloud_path)
        data.text = True
        data.write(output_path + ".txt")
        print("Wrote {}".format(output_path + ".txt"))
    else:
        print("Error:", process.stderr)


def ptcloud_convex_hull_points(pointcloud: np.ndarray) -> np.ndarray:
    assert pointcloud is not None
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=2, tol=1.5, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points


def estimate_normals(convex_hull_surface_pts: np.ndarray, output_path: str, kdtree_radius=None, kdtree_max_nn=None, correction_tangent_planes_n=None): 
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=kdtree_radius, max_nn=kdtree_max_nn) )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    o3d.io.write_point_cloud(output_path, pcd)
    print("Wrote surface with normals {}".format(output_path))
