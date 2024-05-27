import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import ResiduesArg, Chains, ChainArg, StructureArg
from chimerax.atomic import Structure, AtomicStructure, Chain


RIBETL_DATA = os.environ.get("RIBETL_DATA", "/home/rtviii/dev/RIBETL_DATA")


def ribosome_representation(session, structure: AtomicStructure):

    from chimerax.core.commands import run
    from chimerax.core.colors import hex_color

    rcsb_id = str(structure.name).upper()

    with open(os.path.join(RIBETL_DATA, rcsb_id, "{}.json".format(rcsb_id)), "r") as f:
        profile        = json.load(f)


    by_aaid        = {}
    polymer_chains = [ *profile["proteins"], *profile["rnas"], *profile["other_polymers"] ]
    [ by_aaid.update(x) for x in [ {chain['auth_asym_id']:chain} for chain in polymer_chains ] ]

    c:Chain
    for c in structure.chains:
        chain_auth = c.chain_id

        print(chain_auth)
        print(by_aaid[chain_auth])
        # for c in :
        #     nomenclature = c['nomenclature']
        #     if len(nomenclature) < 1:
        #         continue

        #     poly_class = nomenclature[0]
        #     print(poly_class)

    run(session, 'light soft' )
        # run(session, 'color #%s/%s %s ribbon' % (s.id_string, cid, hex_color(chain_colors[cid])))

        #     for r in c["residues"]:
        #         residue = chain.residues[r["residue_id"]]
        #         for a in r["atoms"]:
        #             atom = residue.atoms[a["atom_id"]]
        #             atom.color = a["color"]
        #             atom.display = a["display"]
        #             atom.radius = a["radius"]
        # print(profile)


def register_command(logger):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.core.commands import (
        OpenFolderNameArg,
        BoolArg,
        FloatArg,
        RepeatOf,
        StringArg,
    )
    from chimerax.map import MapArg
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom
    from chimerax.atomic import Residue, Atom

    desc = CmdDesc(
        required=[("structure", AtomicStructureArg)],
        required_arguments=["structure"],
        synopsis="representation ",
    )
    register("ribrep", desc, ribosome_representation, logger=logger)


register_command(session.logger)
