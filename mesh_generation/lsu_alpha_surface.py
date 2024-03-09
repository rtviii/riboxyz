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
from mesh_generation.libsurf import estimate_normals
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, DBSCAN_CLUSTERS_particular_eps_minnbrs
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, normalize_atom_coordinates)
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_chains_by_auth_asym_id
import pyvista as pv
import numpy as np

def lsu_ensemble_get_chains(rcsb_id:str, reconstructed_tunnel_ply:str, outpath:str)->str:
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

def lsu_ensemble_convex_hull(rcsb_id:str, mmcif_ensemble_lsu:str, alpha,tol):
    mmcif_parser     = MMCIFParser(QUIET=True)
    structure        = mmcif_parser.get_structure(rcsb_id+"alpha", mmcif_ensemble_lsu)
    atoms            = np.array([a.get_coord() for a in list(structure.get_atoms())])
    cloud            = pv.PolyData(atoms)
    delaunay_shape   = cloud.delaunay_3d(alpha=alpha, tol=tol, progress_bar=True)
    convex_hull      = delaunay_shape.extract_surface().cast_to_pointset()
    return convex_hull
