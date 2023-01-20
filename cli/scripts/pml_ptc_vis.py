import json
import os
from pprint import pprint
from typing import List
from pymol import cmd

def build_selection_string(site_name: str, chain_name: str, res_ids: List[int]):
    """pymol command being issued"""
    selection_name    = '23S_{}_{}'.format(chain_name, site_name)
    selection_string  = f"c. {chain_name} and "
    selection_string += " OR ".join([*map(lambda x: f"resi {x}", res_ids)])
    return selection_name, selection_string

def highlight_ptc(structure:str):
    structure = structure.upper()
    ptc_coord_path = os.path.join(
                "/home/rxz/dev/docker_ribxz/cli/scripts/PTC_COORDINATES", 
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
        cmd.show("surface", name6)
        cmd.show("surface", name8)
        cmd.show("surface", name9)
        cmd.set("cartoon_transparency", 0.75)
        cmd.reset()

def list_bacteria():
    pprint(os.listdir("/home/rxz/dev/docker_ribxz/cli/scripts/PTC_COORDINATES"))


def ptc(struct:str):
    cmd.delete("all")
    struct      = struct.upper()
    struct_path = os.path.join(
                "/home/rxz/dev/static/{}/{}.cif".format(struct,struct)
            )
    cmd.load(struct_path)
    highlight_ptc(struct)

cmd.extend("ptc",ptc)
cmd.extend("list_bacteria",list_bacteria)