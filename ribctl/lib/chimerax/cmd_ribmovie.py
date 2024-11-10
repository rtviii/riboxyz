import enum
import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript
from chimerax.core.commands import CmdDesc, register, StringArg


def produce_and_save_movie(session, structure:str):
    RCSB_ID = structure
    run(session, "open /home/rtviii/dev/RIBETL_DATA/{}/{}.cif".format(RCSB_ID, RCSB_ID))
    run(session, "sym #1 assembly 1") # take only one assembly if multiple are available
    run(session, "ribrep #2")
    run(session, "movie record")
    run(session, "turn y 2 180")
    run(session, "wait 180")
    run(session, "movie encode /home/rtviii/dev/riboxyz/chimerax/movies/{}.mp4".format(RCSB_ID))
    run(session, "close all")

def register_ribmovie_command(logger):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom

    desc = CmdDesc(
        required           = [("structure", StringArg)],
        required_arguments = ["structure"],
        synopsis           = "making movies",
    )
    register("ribmovie", desc, produce_and_save_movie, logger=logger)


register_ribmovie_command(session.logger)