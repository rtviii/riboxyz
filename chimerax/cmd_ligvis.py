
import enum
import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript
from chimerax.core.commands import CmdDesc, register, StringArg

def register_ligvis_command(logger):
    def ligvis(session, rcsb_id_chemid:str):
        rcsb_id,chemid  = rcsb_id_chemid.split("_")
        profile         = os.path.join("/home/rtviii/dev/riboxyz/antibiotic_bsites", "{}_{}.json".format(rcsb_id,chemid)) 
        run(session, "ribetl {}".format(rcsb_id))
        run(session, "sym #1 assembly 1" )
        run(session, "hide #2" )
        run(session, "show surf")
        run(session, "transparency 85")

        with open(profile, 'r') as f:
            data:dict = json.load(f)

        for auth_asym_id, poly in data.items():
            for resi in poly['residues']:
                run(session, "cartoon #2.1/{}/{} red".format(auth_asym_id,resi['seqid']))
                run(session, "color #2.1/{}/{} red".format(auth_asym_id,resi['seqid']))
        
        run(session, "save /home/rtviii/dev/riboxyz/antibiotic_bsites/images/{}_{}.png".format(rcsb_id,chemid))
        run(session, "close all")

    desc = CmdDesc( required= [("rcsb_id_chemid", StringArg)], required_arguments = ["rcsb_id_chemid"] )
    register("ligvis", desc, ligvis, logger=logger)
