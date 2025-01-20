import os
import sys
from chimerax.core.commands import run ,StringArg
from chimerax.atomic import all_atomic_structures

RIBETL_DATA = os.environ.get("RIBETL_DATA")

def chainsplitter(session, rcsb_id:str):
    run(session, "split chains")
    rcsb_id     = rcsb_id.upper()
    cid_to_spec = {}

    for s in all_atomic_structures(session):
        for chain in s.chains:
            cid_to_spec[chain.chain_id] = s.atomspec

    chains_dir = os.path.join(RIBETL_DATA, rcsb_id, "CHAINS")
    if not os.path.exists(chains_dir):
        os.makedirs(chains_dir)

    for auth_asym_id, model_n in cid_to_spec.items():
        chain_path = os.path.join(chains_dir, f"{rcsb_id}_{auth_asym_id}.cif")
        run(session, "save {} {}".format(chain_path, model_n))
        print("Saved ", chain_path)

def register_ribrepr_command(logger):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom
    desc = CmdDesc(
        required           = [("rcsb_id", StringArg)],
        required_arguments = ["rcsb_id"],
        synopsis           = "representation")
    register("chainsplitter", desc, chainsplitter, logger=logger)

register_ribrepr_command(session.logger)