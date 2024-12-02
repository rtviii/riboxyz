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
from mesh_generation.mesh_bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from mesh_generation.mesh_paths import *
from ribctl import RIBETL_DATA
from ribctl.lib.libpdb import extract_lsu_ensemble_tunnel_vicinity
import pyvista as pv
import numpy as np


"""These methods are for extracting the LSU (or its subset) and clipping the resulting tunnel shape with 
the alphashape of the LSU. I chose to trim manually instead given the complexity of the intersection shapes.
"""

def vestibule_sphere_expansion(rcsb_id:str, radius=50):
    """We want to construct a shape with which to trim the solvent space from obtained tunnel
    Grab the last coordinate in the MOLE centerline and expand to 50 angstrom (captures the vestibule + some interior of the LSU)
    This is an alternative to just constructing a surface from the whole LSU ensemble (LSU rRNA + uL22, uL4, eL39, eL32, uL24) 
    that i'm trying because the guiding arms on the side of the rRNA are tricky to capture with a-shape/normal estimation/poisson recon.(can't obtain a watertight surface)
    """
    [_,x,y,z] = open_tunnel_csv(rcsb_id)[-1]

    mmcif_parser = MMCIFParser(QUIET=True)
    structure    = mmcif_parser.get_structure(rcsb_id, os.path.join(RIBETL_DATA,rcsb_id,rcsb_id + ".cif"))
    atoms        = list(structure.get_atoms())
    ns           = NeighborSearch(atoms)
    _ = ns.search(np.array([x,y,z,]), radius)
    return np.array([a.get_coord() for a in _])

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

    extract_lsu_ensemble_tunnel_vicinity(rcsb_id,list( neighbor_chains_auth_asym_ids ), outpath)
    return outpath

def lsu_ensemble_convex_hull(rcsb_id:str, mmcif_ensemble_lsu:str, alpha,tol):
    mmcif_parser     = MMCIFParser(QUIET=True)
    structure        = mmcif_parser.get_structure(rcsb_id+"alpha", mmcif_ensemble_lsu)
    atoms            = np.array([a.get_coord() for a in list(structure.get_atoms())])
    cloud            = pv.PolyData(atoms)
    delaunay_shape   = cloud.delaunay_3d(alpha=alpha, tol=tol, progress_bar=True)
    convex_hull      = delaunay_shape.extract_surface().cast_to_pointset()
    return convex_hull

def ptcloud_convex_hull( ptcloud:np.ndarray, alpha,tol, offset):
    cloud            = pv.PolyData(ptcloud)
    delaunay_shape   = cloud.delaunay_3d(alpha=alpha, tol=tol, offset=offset,progress_bar=True)
    convex_hull      = delaunay_shape.extract_surface().cast_to_pointset()
    return convex_hull


