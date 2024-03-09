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
from mesh_generation.bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from compas.geometry import bounding_box
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, DBSCAN_CLUSTERS_particular_eps_minnbrs
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, normalize_atom_coordinates)
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_chains_by_auth_asym_id
import pyvista as pv
import numpy as np

def save_lsu_alpha_chains(rcsb_id:str, reconstructed_tunnel_ply:str, outpath:str)->str:
    
    mesh         = pv.read(reconstructed_tunnel_ply)
    points_array = np.array(mesh.points)
    mmcif_parser = MMCIFParser(QUIET=True)
    structure    = mmcif_parser.get_structure(rcsb_id, os.path.join(RIBETL_DATA,rcsb_id,rcsb_id + ".cif"))
    atoms        = list(structure.get_atoms())
    ns           = NeighborSearch(atoms)
    distance_threshold = 5.0
    neighbor_chains_auth_asym_ids = set()
    for point in points_array:
        _ = ns.search(point, distance_threshold)
        [neighbor_chains_auth_asym_ids.add(chain_name) for chain_name in [ a.get_full_id()[2] for a in _]]

    extract_chains_by_auth_asym_id(rcsb_id,list( neighbor_chains_auth_asym_ids ), outpath)
    return outpath


def estimate_normals(rcsb_id, convex_hull_surface_pts: np.ndarray):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=15) )
    pcd.orient_normals_consistent_tangent_plane(k=10)
    o3d.io.write_point_cloud(surface_with_normals_path(rcsb_id), pcd)
    print("Wrote surface with normals {}".format(surface_with_normals_path(rcsb_id)))

def construct_trimming_alphashape(rcsb_id:str, lsu_chains_file:str, alpha,tol):
    mmcif_parser     = MMCIFParser(QUIET=True)
    structure        = mmcif_parser.get_structure(rcsb_id+"alpha", lsu_chains_file)
    atoms            = np.array([a.get_coord() for a in list(structure.get_atoms())])
    cloud            = pv.PolyData(atoms)
    delaunay_shape   = cloud.delaunay_3d(alpha=alpha, tol=tol, progress_bar=True)
    convex_hull      = delaunay_shape.extract_surface().cast_to_pointset()
    estimate_normals(rcsb_id, convex_hull.points)
    return delaunay_shape,convex_hull.points, surface_polydata

def trim_with_alphasurface(rcsb_id:str,reconstruction_pcl:np.ndarray, alpha:float)->np.ndarray:
    """
    Trim the poisson reconstruction with the alpha shape surface
    """
    pass