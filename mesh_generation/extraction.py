
#! ------------------------------ MESH GENERATION
import os
from pprint import pprint
from typing import Tuple
import numpy as np
import open3d as o3d
from open3d import core as o3c
from mendeleev import element
import pandas as pd
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
    ResidueSummary,
)
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from functools import reduce
from Bio.PDB.Structure import Structure
from ribctl.lib.ribosome_types.types_ribosome import RNA


# Tunnel refinement:
# - using the centerline and dynamic probe radius, extract the atoms within 15A radius of the centerline
# - when processing atoms, encode their vdw radius, atom type and residue and chain id
from ribctl import ASSETS_PATH, RIBETL_DATA


def open_tunnel_csv(rcsb_id: str) -> list[list]:
    TUNNEL_PATH = os.path.join(
        ASSETS_PATH, "mole_tunnels", "tunnel_{}.csv".format(rcsb_id)
    )
    df = pd.read_csv(TUNNEL_PATH)
    data = []

    for index, row in df.iterrows():
        radius = row["Radius"]
        x_coordinate = row["X"]
        y_coordinate = row["Y"]
        z_coordinate = row["Z"]
        data.append([radius, x_coordinate, y_coordinate, z_coordinate])

    return data


def parse_struct_via_centerline(rcsb_id: str, centerline_data: list, expansion_radius:int=15) -> list:
    """centerline data is an array of lists [radius, x, y, z]"""
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB import Selection

    parser      = MMCIFParser()
    struct_path = "{}/{}/{}.cif".format(RIBETL_DATA, rcsb_id, rcsb_id)
    structure   = parser.get_structure(rcsb_id, struct_path)
    atoms       = Selection.unfold_entities(structure, "A")
    ns          = NeighborSearch(atoms)
    nbhd        = set()

    for [probe_radius, x, y, z] in centerline_data:
        nearby_atoms = ns.search([x, y, z], probe_radius + expansion_radius, "A")
        nbhd.update(nearby_atoms)

    return list(nbhd)

def parse_struct_via_bbox(rcsb_id: str, bbox:Tuple[list,list]) -> list:
    """bbox is a tuple of minx,miny,minz and maxx,maxy,maxz points"""
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB import Selection

    parser      = MMCIFParser()
    struct_path = "{}/{}/{}.cif".format(RIBETL_DATA, rcsb_id, rcsb_id)
    structure   = parser.get_structure(rcsb_id, struct_path)
    atoms       = Selection.unfold_entities(structure, "A")
    ns          = NeighborSearch(atoms)
    nbhd        = set()

    for [probe_radius, x, y, z] in centerline_data:
        nearby_atoms = ns.search([x, y, z], probe_radius + expansion_radius, "A")
        nbhd.update(nearby_atoms)

    return list(nbhd)

def encode_atoms(rcsb_id: str, nearby_atoms_list: list[Atom]):
    """given a list of atoms lining the tunnel, annotate each with:
    - parent chain id
    - nomenclature
    - containing residue type
    - van der waals radius
    - atom type
    """
    profile      = RibosomeAssets(rcsb_id).profile()
    nomenclature = profile.get_nomenclature_map()
    vdw_radii    = { }
    aggregate = []


    for a in nearby_atoms_list:
        parent_residue      = a.get_parent()
        residue_name        = parent_residue.resname
        residue_seqid       = parent_residue.id[1]
        chain_auth_asym_id  = a.get_full_id()[2]
        parent_nomenclature = nomenclature[chain_auth_asym_id]
        a_element           = a.element
        
        if a_element not in vdw_radii:
            vdw_radii[a_element] = element(a_element).vdw_radius / 100

        atom_dict = { 
                    "coord"             : a.get_coord().tolist(),
                    "chain_auth_asym_id": chain_auth_asym_id,
                    "chain_nomenclature": parent_nomenclature,
                    "residue_name"      : residue_name,
                    "residue_seqid"     : residue_seqid,
                    "atom_element"      : a.element,
                    "vdw_radius"        : vdw_radii[a_element],
                      }
        aggregate.append(atom_dict)
    return aggregate

def create_pcd_from_atoms( positions: np.ndarray, atom_types: np.ndarray, save_path: str ):
    pcd        = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(positions))
    pcd.colors = o3d.utility.Vector3dVector(atom_types)
    o3d.io.write_point_cloud(save_path, pcd)

#TODO : 1. Remove ions and non-standard residues when refining the tunnel
#TODO : CHECK. make the centerline scan in refinement dynamic (on radius of probe) 

def get_sphere_indices_voxelized(center: np.ndarray, radius: int):

    """Make sure radius reflects the size of the underlying voxel grid"""
    x0, y0, z0 = center

    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint
    x_range = slice(int(np.floor(x0) - (radius )), int(np.ceil(x0) + (radius )))
    y_range = slice(int(np.floor(y0) - (radius )), int(np.ceil(y0) + (radius )))
    z_range = slice(int(np.floor(z0) - (radius )), int(np.ceil(z0) + (radius )))

    indices = np.indices((
                x_range.stop - x_range.start,
                y_range.stop - y_range.start,
                z_range.stop - z_range.start
            )
    )

    indices      += np.array([x_range.start, y_range.start, z_range.start])[ :, np.newaxis, np.newaxis, np.newaxis ]
    indices       = indices.transpose(1, 2, 3, 0)
    indices_list  = list(map(tuple, indices.reshape(-1, 3)))
    #!------ Generate indices of a voxel cube of side 2r+2  around the centerpoint

    sphere_active_ix = []

    for ind in indices_list:
        x_ = ind[0]
        y_ = ind[1]
        z_ = ind[2]
        if (x_ - x0) ** 2 + (y_ - y0) ** 2 + (z_ - z0) ** 2 <= radius**2:
            sphere_active_ix.append([x_, y_, z_])

    return np.array(sphere_active_ix)


def centerline_get_bbox(rcsb_id:str):

    data = np.array(open_tunnel_csv(rcsb_id))
    Cx = data[:, 1]
    Cy = data[:, 2]
    Cz = data[:, 3]


    return [ np.min(Cx), np.min(Cy),  np.min(Cz) ],[ np.max(Cx), np.max(Cy), np.max(Cz) ]

print(centerline_get_bbox("6Z6K"))

(bb_p1, bb_p2) = parse_struct_via_centerline("6Z6K", open_tunnel_csv("6Z6K"))