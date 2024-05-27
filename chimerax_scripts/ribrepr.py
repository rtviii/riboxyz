import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import ResiduesArg, Chains, ChainArg, StructureArg
from chimerax.atomic import Structure, AtomicStructure, Chain


RIBETL_DATA = os.environ.get("RIBETL_DATA", "/home/rtviii/dev/RIBETL_DATA")
CHIMERAX_COLORS = [
    ["tan", "#d2b48c"],
    ["sienna", "#a0522d"],
    ["brown", "#a52a2a"],
    ["dark red", "#8b0000"],
    ["firebrick", "#b22222"],
    ["salmon", "#fa8072"],
    ["red", "#ff0000"],
    ["coral", "#ff7f50"],
    ["sandy brown", "#f4a460"],
    ["orange red", "#ff4500"],
    ["orange", "#ff7f00"],
    ["goldenrod", "#daa520"],
    ["gold", "#ffd700"],
    ["yellow", "#ffff00"],
    ["khaki", "#f0e68c"],
    ["dark khaki", "#bdb76b"],
    ["dark olive green", "#556b2f"],
    ["olive drab", "#6b8e23"],
    ["chartreuse", "#7fff00"],
    ["green", "#00ff00"],
    ["dark green", "#006400"],
    ["forest green", "#228b22"],
    ["lime green", "#32cd32"],
    ["light green", "#90ee90"],
    ["sea green", "#2e8b57"],
    ["spring green", "#00ff7f"],
    ["dark cyan", "#008b8b"],
    ["light sea green", "#20b2aa"],
    ["turquoise", "#40e0d0"],
    ["aquamarine", "#7fffd4"],
    ["cyan", "#00ffff"],
    ["deep sky blue", "#00bfff"],
    ["dodger blue", "#1e90ff"],
    ["steel blue", "#4682b4"],
    ["sky blue", "#87ceeb"],
    ["light blue", "#add8e6"],
    ["blue", "#0000ff"],
    ["medium blue", "#3232cd"],
    ["cornflower blue", "#6495ed"],
    ["navy blue", "#000080"],
    ["dark slate blue", "#483d8b"],
    ["medium purple", "#9370db"],
    ["purple", "#a020f0"],
    ["plum", "#dda0dd"],
    ["orchid", "#da70d6"],
    ["magenta", "#ff00ff"],
    ["dark magenta", "#8b008b"],
    ["violet red", "#d02090"],
    ["hot pink", "#ff69b4"],
    ["pink", "#ffc0cb"],
    ["deep pink", "#ff1493"],
    ["rosy brown", "#bc8f8f"],
    ["slate gray", "#708090"],
    ["dark slate gray", "#2f4f4f"],
    ["white", "#ffffff"],
    ["light gray", "#d3d3d3"],
    ["gray", "#bebebe"],
    ["dark gray", "#a9a9a9"],
    ["dim gray", "#696969"],
    ["black", "#000000"],
]


def ribosome_representation(session, structure: AtomicStructure):

    from chimerax.core.commands import run
    from chimerax.core.colors import hex_color
    from chimerax.atomic import Residue, Atom, Chain

    rcsb_id = str(structure.name).upper()

    run(session, "hide #1")

    with open(os.path.join(RIBETL_DATA, rcsb_id, "{}.json".format(rcsb_id)), "r") as f:
        profile = json.load(f)

    polymers = {}
    polymer_chains = [
        *profile["proteins"],
        *profile["rnas"],
        *profile["other_polymers"],
    ]
    [
        polymers.update(x)
        for x in [{chain["auth_asym_id"]: chain} for chain in polymer_chains]
    ]

    c: Chain
    for c in structure.chains:

        aaid = c.chain_id
        polyclass = None

        if len(polymers[aaid]["nomenclature"]) < 1:
            continue
        else:
            polyclass = polymers[aaid]["nomenclature"][0]

        c.polymer_class = polyclass
        Chain.register_attr(session, "polymer_class", "Polymer Class", attr_type=str)
        # Atom.register_attr(session, 'helix_lipophilicity', 'Helix lipophilicity', attr_type = float)

        return

        if polymers[aaid]["entity_poly_polymer_type"] == "RNA":
            run(session, "surf /{}".format(aaid))
            run(session, "color /{} white".format(aaid))
            run(session, "transparency /{} 75".format(aaid))

        run(session, "light soft")
        # for c in :
        #     nomenclature = c['nomenclature']
        #     if len(nomenclature) < 1:
        #         continue

        #     poly_class = nomenclature[0]
        #     print(poly_class)

    run(session, "graphics silhouettes true width 1")
    run(session, "light soft")


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
