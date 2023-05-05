import json
import os
from pprint import pprint
from typing import List
from api.scratch_tunnel_workflow import ptc_midpoint, ptc_residues_via_alignment, tunnel_obstructions
from pymol import cmd

RIBETL_DATA = os.environ.get("RIBETL_DATA")

def build_selection_string( chain_name: str, res_ids: List[int]):
    """pymol command being issued"""
    selection_name    = '23-28S_{}'.format(chain_name)
    selection_string  = f"c. {chain_name} and "
    selection_string += " OR ".join([*map(lambda x: f"resi {x}", res_ids)])
    return selection_name, selection_string

def list_bacteria():
    pprint(os.listdir(f"{RIBETL_DATA}/PTC_COORDINATES"))

def ptc(struct: str):
    cmd.delete("all")
    struct = struct.upper()
    struct_path = os.path.join("/home/rxz/dev/static/{}/{}.cif".format(struct, struct))
    cmd.load(struct_path)
    highlight_ptc(struct)

def get_markerspath(struct: str):
    struct = struct.upper()
    _path = os.path.join(RIBETL_DATA, "PTC_MARKERS_RAW", f"{struct}_PTC_MARKERS_RAW.json")
    return _path

def create_marker_at_atom(selection_name:str, posn:List[float], color_:str="red", repr="spheres", label=''):
    cmd.pseudoatom(selection_name, pos=posn, vdw=1, color=color_,  label=label)
    cmd.show(repr, selection_name)

def ptc_raw_w_markerks(struct: str):
    # cmd.delete("all")
    # struct      = struct.upper()
    # struct_path = os.path.join("/home/rxz/dev/static/{}/{}.cif".format(struct, struct))
    # cmd.load(struct_path)

    highlight_ptc_raw(struct)

    with open(get_markerspath(struct), 'r') as infile:
       POSNS:dict = json.load(infile)

    chain = [*POSNS.keys()][0]
    if chain == None:
        exit("Could not identify chain")

    U_end_pos   = [ *POSNS[chain].values() ][len(POSNS[chain]) - 2]["O4'"] # pre  last residue of the comb
    U_start_pos = [ *POSNS[chain].values() ][0]["O4'"] # firs      residue of the comb

    midpoint = [
        ( U_end_pos[0] + U_start_pos[0] ) / 2,
        ( U_end_pos[1] + U_start_pos[1] ) / 2,
        ( U_end_pos[2] + U_start_pos[2] ) / 2,
    ]
    create_marker_at_atom("uridine_comb_start", U_start_pos, color_="green")
    create_marker_at_atom("uridine_comb_end", U_end_pos, color_="green")
    cmd.distance(None, "uridine_comb_start", "uridine_comb_end", mode=0)

    create_marker_at_atom("centroid",midpoint, color_="red")

def visualize_obstructions(rcsb_id):
    ptcres, auth_asym_id = ptc_residues_via_alignment(rcsb_id, 0)
    midpoint = ptc_midpoint(ptcres, auth_asym_id)
    polys,nonpolys = tunnel_obstructions(rcsb_id, midpoint)

    cmd.color("white", "all")

    ptc_raw_w_markerks(rcsb_id)

    for poly in polys:
        cname = "chain_{}".format(poly.auth_asym_id, )
        csele = "c. {} and m. {}".format(poly.auth_asym_id, rcsb_id.upper()) 
        cmd.create(cname, csele)
        cmd.remove(csele)
        cmd.color("red", cname)

    for res in nonpolys:
        resname = "res_{}_{}".format(res.parent_auth_asym_id,res.resname)
        ressele = "m. {} and c. {} and resn {}".format(rcsb_id.upper(),res.parent_auth_asym_id,res.resname)
        cmd.create(resname,  ressele)
        cmd.remove(ressele)
        cmd.color("blue", resname)

# def struct_paint_chains(rcsb_id):

def struct_paint_chains(pdbid: str):
    pdbid          = pdbid.upper()
    RIBETL_DATA    = str(os.environ.get('RIBETL_DATA'))
    nomenclaturev2 = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}_nomenclaturev2.json")
    profilepath    = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.json")

    with open(profilepath, 'rb') as infile:
        profile = json.load(infile)

    if profile['rnas'] != None:
        for rna in profile['rnas']:
            # cmd.color('white', f"chain {rna['auth_asym_id']}")
            cmd.hide('everything', f"chain {rna['auth_asym_id']}")
            cmd.show('surface', f"chain {rna['auth_asym_id']}")
            cmd.color('white', f"chain {rna['auth_asym_id']}")
            cmd.set('transparency', 0.1, f"chain {rna['auth_asym_id']}")

    for protein in profile['proteins']:

        if len(protein['nomenclature']) != 0:
            if protein['nomenclature'][0] in colormap__LSU_Proteins:
                CLR = colormap__LSU_Proteins[protein['nomenclature'][0]]
            elif protein['nomenclature'][0] in colormap__SSU_Proteins:
                CLR = colormap__SSU_Proteins[protein['nomenclature'][0]]
        else:
            CLR = 'blue'

        prot_tmp = f"protein_tmp_{ protein['auth_asym_id'] }"

        cmd.create(prot_tmp, f"chain {protein[ 'auth_asym_id' ]}")
        cmd.remove(f"chain {protein[ 'auth_asym_id' ]} and m. {pdbid}")
        cmd.hide ('everything', prot_tmp)
        cmd.show ('surface', prot_tmp)
        cmd.show ('sticks', prot_tmp)
        cmd.show ('cartoon', prot_tmp)
        cmd.show ('lines', prot_tmp)
        cmd.set('transparency', 0.75, prot_tmp)
        cmd.color(CLR , prot_tmp)

def sload(pdbid: str):
    pdbid       = pdbid.upper()
    RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
    path        = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif")
    cmd.delete('all')
    cmd.load(path)

cmd.extend("sload", sload)
cmd.extend("by_rna", struct_paint_chains)
cmd.extend("ptc", ptc)
cmd.extend("ptc_w_markers", ptc_raw_w_markerks)
cmd.extend("list_bacteria", list_bacteria)
cmd.extend("tun_obstructions", visualize_obstructions)


