#! ------------------------------ MESH GENERATION
import pyvista as pv
import json
import os
from pprint import pprint
import sys
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
from compas.geometry import bounding_box, oriented_bounding_box_numpy

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


def parse_struct_via_centerline(
    rcsb_id: str, centerline_data: list, expansion_radius: int = 15
) -> list[Atom]:
    """centerline data is an array of lists [radius, x, y, z]"""
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB import Selection

    parser = MMCIFParser()
    struct_path = "{}/{}/{}.cif".format(RIBETL_DATA, rcsb_id, rcsb_id)
    structure = parser.get_structure(rcsb_id, struct_path)
    atoms = Selection.unfold_entities(structure, "A")
    ns = NeighborSearch(atoms)
    nbhd = set()

    for [probe_radius, x, y, z] in centerline_data:
        nearby_atoms = ns.search([x, y, z], probe_radius + expansion_radius, "A")
        nbhd.update(nearby_atoms)

    return list(nbhd)


def parse_struct_via_bbox(rcsb_id: str, bbox: list) -> list:
    """bbox is a tuple of minx,miny,minz and maxx,maxy,maxz points"""
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.NeighborSearch import NeighborSearch
    from Bio.PDB import Selection


    parser      = MMCIFParser()
    struct_path = "{}/{}/{}.cif".format(RIBETL_DATA, rcsb_id, rcsb_id)
    structure   = parser.get_structure(rcsb_id, struct_path)
    atoms       = Selection.unfold_entities(structure, "A")
    nbhd        = []

    def is_inside_box(point, box_coordinates):
        EXPANSION_RADIUS = 0

        min_x = min(box_coordinates, key=lambda x: x[0])[0]
        max_x = max(box_coordinates, key=lambda x: x[0])[0]

        min_y = min(box_coordinates, key=lambda x: x[1])[1]
        max_y = max(box_coordinates, key=lambda x: x[1])[1]

        min_z = min(box_coordinates, key=lambda x: x[2])[2]
        max_z = max(box_coordinates, key=lambda x: x[2])[2]

        if  ( min_x+EXPANSION_RADIUS ) <= point[0] <= ( max_x+EXPANSION_RADIUS ) and \
            ( min_y+EXPANSION_RADIUS ) <= point[1] <= ( max_y+EXPANSION_RADIUS ) and \
            ( min_z+EXPANSION_RADIUS ) <= point[2] <= ( max_z+EXPANSION_RADIUS ) :
            return True
        else:
            return False

    for atom in atoms:
        if is_inside_box(atom.get_coord(), bbox):
            nbhd.append(atom)
        

    return list(nbhd)


def encode_atoms(rcsb_id: str, nearby_atoms_list: list[Atom], write=False) -> list:
    """given a list of atoms lining the tunnel, annotate each with:
    - parent chain id
    - nomenclature
    - containing residue type
    - van der waals radius
    - atom type
    """
    profile      = RibosomeAssets(rcsb_id).profile()
    nomenclature = profile.get_nomenclature_map()
    vdw_radii    = {}
    aggregate    = []

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
            "coord": a.get_coord().tolist(),
            "chain_auth_asym_id": chain_auth_asym_id,
            "chain_nomenclature": parent_nomenclature,
            "residue_name": residue_name,
            "residue_seqid": residue_seqid,
            "atom_element": a.element,
            "vdw_radius": vdw_radii[a_element],
        }
        aggregate.append(atom_dict)

    if write:
        with open( "/home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms_bbox.json".format( rcsb_id ), "w", ) as outfile:
            json.dump(aggregate, outfile, indent=4)
        print( "Wrote {} tunnel atoms to disk : /home/rtviii/dev/riboxyz/mesh_generation/{}_tunnel_atoms_bbox.json".format( len(aggregate), rcsb_id ) )

    return aggregate

def create_pcd_from_atoms(
    positions: np.ndarray, atom_types: np.ndarray, save_path: str
):
    pcd = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(positions))
    pcd.colors = o3d.utility.Vector3dVector(atom_types)
    o3d.io.write_point_cloud(save_path, pcd)



RCSB_ID      = sys.argv[1].upper()
IF_VISUALIZE = sys.argv[2].upper()


cloud            = parse_struct_via_centerline(RCSB_ID, open_tunnel_csv(RCSB_ID))
cords_walls_only = np.array([a.get_coord() for a in cloud])

if IF_VISUALIZE:
    point_cloud           = pv.PolyData(cords_walls_only)
    # Coloration
    random_rgbs           = np.random.randint(0, 256, size=( cords_walls_only.shape[0],4 ))
    point_cloud['colors'] = random_rgbs
    point_cloud.plot(scalars='colors', rgb=True, notebook=False)

bbox  = bounding_box(np.array([a.get_coord() for a in cloud]))
print("Vanilla bounding box:", bbox)
atoms_bboxed = parse_struct_via_bbox(RCSB_ID, bbox)

if IF_VISUALIZE:
    cords_inside_bbox = np.array([a.get_coord() for a in atoms_bboxed])
    point_cloud           = pv.PolyData(cords_inside_bbox)
    random_rgbs           = np.random.randint(0, 256, size=( cords_inside_bbox.shape[0],4 ))
    point_cloud['colors'] = random_rgbs
    point_cloud.plot(scalars='colors', rgb=True, notebook=False)

encode_atoms(RCSB_ID, atoms_bboxed, write=True)
