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
from ribctl import ASSETS_PATH, EXIT_TUNNEL_WORK


RCSB_ID = "6Z6K"

tunnel_atom_encoding_handle = "{}_tunnel_atoms_bbox.json".format(RCSB_ID)
tunnel_atom_encoding_path   = os.path.join(EXIT_TUNNEL_WORK,tunnel_atom_encoding_handle)




centerline_expansion_atoms       = parse_struct_via_centerline(RCSB_ID, open_tunnel_csv(RCSB_ID))
centerline_expansion_coordinates = np.array([a.get_coord() for a in centerline_expansion_atoms])
bbox                             = bounding_box(centerline_expansion_coordinates)
bbox_interior_atoms              = parse_struct_via_bbox(RCSB_ID, bbox)

encode_atoms(RCSB_ID, bbox_interior_atoms, write=True)