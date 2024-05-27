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

    run(session, 'hide #1')

    with open(os.path.join(RIBETL_DATA, rcsb_id, "{}.json".format(rcsb_id)), "r") as f:
        profile        = json.load(f)


    polymers        = {}
    polymer_chains = [ *profile["proteins"], *profile["rnas"], *profile["other_polymers"] ]
    [ polymers.update(x) for x in [ {chain['auth_asym_id']:chain} for chain in polymer_chains ] ]

    c:Chain
    for c in structure.chains:
        aaid = c.chain_id

        if polymers[aaid]["entity_poly_polymer_type"] == "RNA":
            run(session, 'surf /{}'.format(aaid))
            run(session, 'color /{} white'.format(aaid))
            run(session, 'transparency /{} 75'.format(aaid))



        run(session, 'light soft' )
        # for c in :
        #     nomenclature = c['nomenclature']
        #     if len(nomenclature) < 1:
        #         continue

        #     poly_class = nomenclature[0]
        #     print(poly_class)

    run(session, 'graphics silhouettes true width 1')
    run(session, 'light soft' )


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
