#! ------------------------------ MESH GENERATION
import json
import os
from pprint import pprint
import numpy as np
import open3d as o3d
from mendeleev import element
import pandas as pd
from ribctl.etl.etl_ribosome_ops import RibosomeOps, Structure
from Bio.PDB.Atom import Atom

# Tunnel refinement:
# - using the centerline and dynamic probe radius, extract the atoms within 15A radius of the centerline
# - when processing atoms, encode their vdw radius, atom type and residue and chain id
from ribctl import  EXIT_TUNNEL_WORK, RIBETL_DATA

def open_tunnel_csv(rcsb_id: str) -> list[list]:
    TUNNEL_PATH = os.path.join( EXIT_TUNNEL_WORK, "mole_tunnels", "tunnel_{}.csv".format(rcsb_id) )
    df          = pd.read_csv(TUNNEL_PATH)
    data        = []
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

    parser      = MMCIFParser()
    struct_path = RibosomeOps(rcsb_id).paths.cif
    structure   = parser.get_structure(rcsb_id, struct_path)
    atoms       = Selection.unfold_entities(structure, "A")
    ns          = NeighborSearch(atoms)
    nbhd        = set()

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
        
        [
            [min_x, min_y, min_z],
            [max_x, min_y, min_z],
            [max_x, max_y, min_z],
            [min_x, max_y, min_z],
            [min_x, min_y, max_z],
            [max_x, min_y, max_z],
            [max_x, max_y, max_z],
            [min_x, max_y, max_z],
                                    ] = box_coordinates


        if  ( min_x ) <= point[0] <= ( max_x ) and \
            ( min_y ) <= point[1] <= ( max_y ) and \
            ( min_z ) <= point[2] <= ( max_z ) :
            return True
        else:
            return False

    for atom in atoms:
        if is_inside_box(atom.get_coord(), bbox):
            nbhd.append(atom)
    return list(nbhd)

def encode_atoms(rcsb_id: str, atoms_list: list[Atom], write=False, writepath=None) -> list:
    """given a list of atoms lining the tunnel, annotate each with:
    - parent chain id
    - nomenclature
    - containing residue type
    - van der waals radius
    - atom type
    """
    profile      = RibosomeOps(rcsb_id).profile()
    nomenclature = profile.get_nomenclature_map()
    vdw_radii    = {}
    aggregate    = []

    for a in atoms_list:
        parent_residue      = a.get_parent()
        residue_name        = parent_residue.resname
        residue_seqid       = parent_residue.id[1]
        chain_auth_asym_id  = a.get_full_id()[2]
        parent_nomenclature = nomenclature[chain_auth_asym_id]
        a_element           = a.element

        try: 
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
        except Exception as e:
            print(f"Couldn't figure out atom {a} :", e)
            print("Skipping...")


    if write and writepath:
        with open( writepath, "w", ) as outfile:
            json.dump(aggregate, outfile, indent=4)
        print("Wrote {} tunnel atoms to disk at {}".format( len(aggregate), writepath ) )
    
    if write and not writepath:
        raise LookupError("Provide writepath to `encode_atoms`.")

    return aggregate

def create_pcd_from_atoms( positions: np.ndarray, atom_types: np.ndarray, save_path: str ):
    pcd        = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(positions))
    pcd.colors = o3d.utility.Vector3dVector(atom_types)
    o3d.io.write_point_cloud(save_path, pcd)

