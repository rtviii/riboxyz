import json
import os
from chimerax.core.commands import register, CmdDesc
from chimerax.core.commands import run
from chimerax.core.commands import CmdDesc, register, StringArg

def register_ligvis_command(logger):
    def ligvis(session, rcsb_id_chemid:str):
        chemid, rcsb_id  = rcsb_id_chemid.split("_")
        OUTPATH = "/home/rtviii/dev/riboxyz/antibiotic_bsites/images/{}_{}.png".format(chemid, rcsb_id)
        if os.path.exists(OUTPATH):
            return

        run(session, "set bgColor white")
        profile         = os.path.join("/home/rtviii/dev/riboxyz/antibiotic_bsites", "{}_{}.json".format(chemid, rcsb_id)) 
        run(session, "ribetl {}".format(rcsb_id))
        run(session, "sym #1 assembly 1" )
        run(session, "hide #2" )
        run(session, "show surf")
        run(session, "transparency 85")

        with open(profile, 'r') as f:
            data:dict = json.load(f)

        for auth_asym_id, poly in data.items():
            for resi in poly['residues']:
                run(session, "cartoon #2/{}:{}".format(auth_asym_id,resi['seqid']))
                run(session, "color #2/{}:{} red".format(auth_asym_id,resi['seqid']))
        
        run(session, "save /home/rtviii/dev/riboxyz/antibiotic_bsites/images/{}_{}.png width 800 height 800".format(chemid, rcsb_id))
        run(session, "close all")
        print("Produced image for ", chemid, rcsb_id)

    desc = CmdDesc( required= [("rcsb_id_chemid", StringArg)], required_arguments = ["rcsb_id_chemid"] )
    register("ligvis", desc, ligvis, logger=logger)


register_ligvis_command(session.logger)