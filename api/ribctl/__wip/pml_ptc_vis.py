import json
import os
from pprint import pprint
from typing import List
from pymol import cmd


RIBETL_DATA = "/home/rxz/dev/static"


def build_selection_string( chain_name: str, res_ids: List[int]):
    """pymol command being issued"""
    selection_name    = '23-28S_{}'.format(chain_name)
    selection_string  = f"c. {chain_name} and "
    selection_string += " OR ".join([*map(lambda x: f"resi {x}", res_ids)])
    return selection_name, selection_string


def highlight_ptc(structure: str):
    structure = structure.upper()
    ptc_coord_path = os.path.join(
        RIBETL_DATA,
        "PTC_COORDINATES",
        f'{structure}_PTC.json'
    )
    with open(ptc_coord_path, 'r') as f:
        ptc = json.load(f)
        chain = [*ptc.keys()][0]
        if chain == None:
            exit("Could not identify chain")

        for site in ["site_6", "site_8", "site_9"]:
            name6, selection6 = build_selection_string(
                site,
                chain,
                [*map(lambda truple: truple[0], ptc[chain][site])]
            )

        name6, selection6 = build_selection_string(
            "site_6", chain, [*map(lambda truple: truple[0], ptc[chain]["site_6"])])
        name8, selection8 = build_selection_string(
            "site_8", chain, [*map(lambda truple: truple[0], ptc[chain]["site_8"])])
        name9, selection9 = build_selection_string(
            "site_9", chain, [*map(lambda truple: truple[0], ptc[chain]["site_9"])])

        cmd.create(name6, selection6)
        cmd.create(name8, selection8)
        cmd.create(name9, selection9)

        cmd.color("gray60", structure)

        cmd.color("forest", name6)
        cmd.color("tv_orange", name8)
        cmd.color("blue", name9)

        cmd.show("surface", name6)
        cmd.show("surface", name8)
        cmd.show("surface", name9)
        # cmd.set("cartoon_transparency", 0.5)
        cmd.bg_color("white")
        cmd.reset()


def highlight_ptc_fuzzy(structure: str):
    """site 9 conserved residues located via fuzzy search"""
    structure = structure.upper()
    ptc_fuzzy_coord_path = os.path.join(
        RIBETL_DATA,
        "PTC_COORDINATES",
        f'{structure}_FUZZY_PTC.json'
    )
    print(ptc_fuzzy_coord_path)
    with open(ptc_fuzzy_coord_path, 'r') as f:
        ptc = json.load(f)
        chain = [*ptc.keys()][0]
        if chain == None:
            exit("Could not identify chain")

        # name6, selection6 = build_selection_string(
        #     "site_6", chain, ptc[chain]["site_6"])
        # name8, selection8 = build_selection_string(
        #     "site_8", chain, ptc[chain]["site_8"])
        name9, selection9 = build_selection_string("site_9", chain, ptc[chain]["site_9"])

        # cmd.create(name6, selection6)
        # cmd.create(name8, selection8)
        cmd.create(name9, selection9)

        cmd.color("gray60", structure)

        # cmd.color("forest", name6)
        # cmd.color("tv_orange", name8)

        # cmd.show("surface", name6)
        # cmd.show("surface", name8)
        # cmd.show("surface", name9)
        # cmd.show("licorice", name6)
        # cmd.show("licorice", name8)
        # cmd.set("cartoon_transparency", 0.5)
        cmd.bg_color("white")
        cmd.set("transparency", 0.5)


        # SITE_9_OBJ ="site9"
        # cmd.create(SITE_9_OBJ, name9)

        cmd.show("licorice", name9)
        cmd.color("blue", name9)
        # cmd.label(name9, "resi")
        cmd.set("label_color", "pink")
        cmd.set("label_size", 5)


        cmd.zoom(name9)
        cmd.reset()

def highlight_ptc_raw(structure: str):
    structure            = structure.upper()
    ptc_fuzzy_coord_path = os.path.join(
        RIBETL_DATA,
        "PTC_MARKERS_RAW",
        f'{structure}_PTC_MARKERS_RAW.json'
    )
    print(ptc_fuzzy_coord_path)
    with open(ptc_fuzzy_coord_path, 'r') as f:
        ptc = json.load(f)
        chain = [*ptc.keys()][0]
        if chain == None:
            exit("Could not identify chain")

        name9, selection9 = build_selection_string(chain, ptc[chain].keys())

        cmd.create(name9, selection9)
        cmd.color("gray60", structure)

        # cmd.set("cartoon_transparency", 0.5)
        cmd.bg_color("white")
        cmd.set("transparency", 0.5)

        cmd.show("licorice", name9)
        cmd.color("blue", name9)
        cmd.set("label_color", "pink")
        cmd.set("label_size", 5)

        cmd.zoom(name9)
        cmd.reset()


def list_bacteria():
    pprint(os.listdir(f"{RIBETL_DATA}/PTC_COORDINATES"))


def ptc(struct: str):
    cmd.delete("all")
    struct = struct.upper()
    struct_path = os.path.join(
        "/home/rxz/dev/static/{}/{}.cif".format(struct, struct))
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
    cmd.delete("all")
    struct      = struct.upper()
    struct_path = os.path.join("/home/rxz/dev/static/{}/{}.cif".format(struct, struct))
    cmd.load(struct_path)

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
    # print( "PRESENT::",[ *POSNS[chain].values() ][0].keys())
    create_marker_at_atom("uridine_comb_start",U_start_pos, color_="green")
    create_marker_at_atom("uridine_comb_end",U_end_pos, color_="green")
    cmd.distance(None, "uridine_comb_start", "uridine_comb_end", mode=0)

    create_marker_at_atom("centroid",midpoint, color_="red")

cmd.extend("ptc", ptc_raw_w_markerks)
cmd.extend("list_bacteria", list_bacteria)


