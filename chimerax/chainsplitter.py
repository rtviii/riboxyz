import enum
import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import ResiduesArg, Chains, ChainArg, StructureArg
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript
from chimerax.list_info.util import spec, model_info
from chimerax.core.state import StateManager

def chainsplitter(session, structure: AtomicStructure):
    print("Got some chains")
    print("----------------------------------")
    # print(report_polymers(session.logger, structure))
    # run(session, "split chains")
    # chain:Chain
    print(model_info(structure))
    # for model in all
    #     print(model)
    # for chain in  structure.chains:
    #     run(session,"save ")
    #     print(chain.chain_id)
    # print(model_info(structure))
    

def register_ribrepr_command(logger):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom

    desc = CmdDesc(
        required           = [("structure", AtomicStructureArg)],
        required_arguments = ["structure"],
        synopsis           = "representation ",
    )
    register("chainsplitter", desc, chainsplitter, logger=logger)

register_ribrepr_command(session.logger)