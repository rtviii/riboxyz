import enum
import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript
from chimerax.core.commands import CmdDesc, register, StringArg


RIBETL_DATA = os.environ.get("RIBETL_DATA")

def register_ribetl_command(logger):
    def ribetl(session, rcsb_id:str):
        rcsb_id = rcsb_id.upper()
        cifpath = os.path.join(RIBETL_DATA, rcsb_id, "{}.cif".format(rcsb_id))
        run(session, "open {}".format(rcsb_id,cifpath))
    desc = CmdDesc( required= [("rcsb_id", StringArg)], required_arguments = ["rcsb_id"] )
    register("ribetl", desc, ribetl, logger=logger)

register_ribetl_command(session.logger)