from pprint import pprint
import open3d as o3d
import json
import os
import numpy as np
from mesh_generation.bbox_extraction import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_bbox,
    parse_struct_via_centerline,
)
import concurrent.futures

from compas.geometry import bounding_box, oriented_bounding_box_numpy
from mesh_generation.voxelize import expand_atomcenters_to_spheres_threadpool, get_sphere_indices_voxelized
from ribctl import ASSETS_PATH, EXIT_TUNNEL_WORK


RCSB_ID = "6Z6K"

tunnel_atom_encoding_handle = "{}_tunnel_atoms_bbox.json".format(RCSB_ID)
tunnel_atom_encoding_path   = os.path.join(EXIT_TUNNEL_WORK,tunnel_atom_encoding_handle)


def extract_bbox_atoms(rcsb_id:str):

    centerline_expansion_atoms       = parse_struct_via_centerline(rcsb_id, open_tunnel_csv(rcsb_id))
    centerline_expansion_coordinates = np.array([a.get_coord() for a in centerline_expansion_atoms])
    bbox                             = bounding_box(centerline_expansion_coordinates)
    bbox_interior_atoms              = parse_struct_via_bbox(rcsb_id, bbox)

    encode_atoms(rcsb_id, bbox_interior_atoms, write=True, writepath=tunnel_atom_encoding_path)


def expand_bbox_atoms_to_spheres(rcsb_id:str):

    with open( tunnel_atom_encoding_path, "r", ) as infile:
        bbox_data = json.load(infile)

    __cords = np.array(list(map(lambda x: x["coord"], bbox_data)))
    __radii = np.array(list(map(lambda x: x["vdw_radius"], bbox_data)))

    sphere_sources = zip(__cords, __radii)

    SINK     = []
    expanded = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)

    return np.array(expanded)