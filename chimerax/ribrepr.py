import enum
import json
import os
import sys
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import ResiduesArg, Chains, ChainArg, StructureArg
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript

# https://mail.cgl.ucsf.edu/mailman/archives/list/chimera-users@cgl.ucsf.edu/thread/EOUA5K3CZU6DJYVPISR3GFWHCISR6WGV/


RIBETL_DATA = os.environ.get("RIBETL_DATA", "/home/rtviii/dev/RIBETL_DATA")

class PolymerClass(str, enum.Enum):

    tRNA = "tRNA"
    bS1m  = "bS1m"
    uS2m  = "uS2m"
    uS3m  = "uS3m"
    uS4m  = "uS4m"
    uS5m  = "uS5m"
    bS6m  = "bS6m"
    uS7m  = "uS7m"
    uS8m  = "uS8m"
    uS9m  = "uS9m"
    uS10m = "uS10m"
    uS11m = "uS11m"
    uS12m = "uS12m"
    uS13m = "uS13m"
    uS14m = "uS14m"
    uS15m = "uS15m"
    bS16m = "bS16m"
    uS17m = "uS17m"
    bS18m = "bS18m"
    uS19m = "uS19m"
    bS21m = "bS21m"
    mS22  = "mS22"
    mS23  = "mS23"
    mS25  = "mS25"
    mS26  = "mS26"
    mS27  = "mS27"
    mS29  = "mS29"
    mS31  = "mS31"
    mS33  = "mS33"
    mS34  = "mS34"
    mS35  = "mS35"
    mS37  = "mS37"
    mS38  = "mS38"
    mS39  = "mS39"
    mS40  = "mS40"
    mS41  = "mS41"
    mS42  = "mS42"
    mS43  = "mS43"
    mS44  = "mS44"
    mS45  = "mS45"
    mS46  = "mS46"
    mS47  = "mS47"
    uL1m  = "uL1m"
    uL2m  = "uL2m"
    uL3m  = "uL3m"
    uL4m  = "uL4m"
    uL5m  = "uL5m"
    uL6m  = "uL6m"
    bL9m  = "bL9m"
    uL10m = "uL10m"
    uL11m = "uL11m"
    bL12m = "bL12m"
    uL13m = "uL13m"
    uL14m = "uL14m"
    uL15m = "uL15m"
    uL16m = "uL16m"
    bL17m = "bL17m"
    uL18m = "uL18m"
    bL19m = "bL19m"
    bL20m = "bL20m"
    bL21m = "bL21m"
    uL22m = "uL22m"
    uL23m = "uL23m"
    uL24m = "uL24m"
    bL27m = "bL27m"
    bL28m = "bL28m"
    uL29m = "uL29m"
    uL30m = "uL30m"
    bL31m = "bL31m"
    bL32m = "bL32m"
    bL33m = "bL33m"
    bL34m = "bL34m"
    bL35m = "bL35m"
    bL36m = "bL36m"
    mL37  = "mL37"
    mL38  = "mL38"
    mL39  = "mL39"
    mL40  = "mL40"
    mL41  = "mL41"
    mL42  = "mL42"
    mL43  = "mL43"
    mL44  = "mL44"
    mL45  = "mL45"
    mL46  = "mL46"
    mL48  = "mL48"
    mL49  = "mL49"
    mL50  = "mL50"
    mL51  = "mL51"
    mL52  = "mL52"
    mL53  = "mL53"
    mL54  = "mL54"
    mL57  = "mL57"
    mL58  = "mL58"
    mL59  = "mL59"
    mL60  = "mL60"
    mL61  = "mL61"
    mL62  = "mL62"
    mL63  = "mL63"
    mL64  = "mL64"
    mL65  = "mL65"
    mL66  = "mL66"
    mL67  = "mL67"
    bS1   = "bS1"
    eS1   = "eS1"
    uS2   = "uS2"
    uS3   = "uS3"
    uS4   = "uS4"
    eS4   = "eS4"
    uS5   = "uS5"
    bS6   = "bS6"
    eS6   = "eS6"
    uS7   = "uS7"
    eS7   = "eS7"
    uS8   = "uS8"
    eS8   = "eS8"
    uS9   = "uS9"
    uS10  = "uS10"
    eS10  = "eS10"
    uS11  = "uS11"
    uS12  = "uS12"
    eS12  = "eS12"
    uS13  = "uS13"
    uS14  = "uS14"
    uS15  = "uS15"
    bS16  = "bS16"
    uS17  = "uS17"
    eS17  = "eS17"
    bS18  = "bS18"
    uS19  = "uS19"
    eS19  = "eS19"
    bS20  = "bS20"
    bS21  = "bS21"
    bTHX  = "bTHX"
    eS21  = "eS21"
    eS24  = "eS24"
    eS25  = "eS25"
    eS26  = "eS26"
    eS27  = "eS27"
    eS28  = "eS28"
    eS30  = "eS30"
    eS31  = "eS31"
    RACK1 = "RACK1"
    uL1  = "uL1"
    uL2  = "uL2"
    uL3  = "uL3"
    uL4  = "uL4"
    uL5  = "uL5"
    uL6  = "uL6"
    eL6  = "eL6"
    eL8  = "eL8"
    bL9  = "bL9"
    uL10 = "uL10"
    uL11 = "uL11"
    bL12 = "bL12"
    uL13 = "uL13"
    eL13 = "eL13"
    uL14 = "uL14"
    eL14 = "eL14"
    uL15 = "uL15"
    eL15 = "eL15"
    uL16 = "uL16"
    bL17 = "bL17"
    uL18 = "uL18"
    eL18 = "eL18"
    bL19 = "bL19"
    eL19 = "eL19"
    bL20 = "bL20"
    eL20 = "eL20"
    bL21 = "bL21"
    eL21 = "eL21"
    uL22 = "uL22"
    eL22 = "eL22"
    uL23 = "uL23"
    uL24 = "uL24"
    eL24 = "eL24"
    bL25 = "bL25"
    bL27 = "bL27"
    eL27 = "eL27"
    bL28 = "bL28"
    eL28 = "eL28"
    uL29 = "uL29"
    eL29 = "eL29"
    uL30 = "uL30"
    eL30 = "eL30"
    bL31 = "bL31"
    eL31 = "eL31"
    bL32 = "bL32"
    eL32 = "eL32"
    bL33 = "bL33"
    eL33 = "eL33"
    bL34 = "bL34"
    eL34 = "eL34"
    bL35 = "bL35"
    bL36 = "bL36"
    eL36 = "eL36"
    eL37 = "eL37"
    eL38 = "eL38"
    eL39 = "eL39"
    eL40 = "eL40"
    eL41 = "eL41"
    eL42 = "eL42"
    eL43 = "eL43"
    P1P2 = "P1P2"
    mtrRNA12S = "mt12SrRNA"  # mitochondrial
    mtrRNA16S = "mt16SrRNA"  # mitochondrial
    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA"  #  c-bacterial or mitochondrial
    rRNA_23S  = "23SrRNA"  # bacterial
    rRNA_25S  = "25SrRNA"  # plants
    rRNA_5_8S = "5.8SrRNA"  # eukaryotic
    rRNA_18S  = "18SrRNA"  # eukaryotic
    rRNA_28S  = "28SrRNA"  # eukaryotic
    # Eukaryotic
    eEF1A = "eEF1A"
    eEF1B = "eEF1B"
    eFSec = "eFSec"
    eEF2  = "eEF2"
    mtEF4 = "mtEF4"
    eIF5A = "eIF5A"
    eEF3  = "eEF3"
    # Bacterial
    EF_Tu = "EF-Tu"
    EF_Ts = "EF-Ts"
    SelB  = "SelB"
    EF_G  = "EF-G"
    EF4   = "EF4"
    EF_P  = "EF-P"
    Tet_O = "Tet_O"
    Tet_M = "Tet_M"
    RelA  = "RelA"
    BipA  = "BipA"
    # Archaeal
    aEF1A = "aEF1A"
    aEF2  = "aEF2"
    #!Eukaryotic
    eIF1  = "eIF1"
    eIF1A = "eIF1A"

    eIF2_alpha = "eIF2_alpha"
    eIF2_beta  = "eIF2_beta"
    eIF2_gamma = "eIF2_gamma"

    eIF2B_alpha   = "eIF2B_alpha"
    eIF2B_beta    = "eIF2B_beta"
    eIF2B_gamma   = "eIF2B_gamma"
    eIF2B_delta   = "eIF2B_delta"
    eIF2B_epsilon = "eIF2B_epsilon"

    eIF3_subunitA = "eIF3_subunitA"
    eIF3_subunitB = "eIF3_subunitB"
    eIF3_subunitC = "eIF3_subunitC"
    eIF3_subunitD = "eIF3_subunitD"
    eIF3_subunitE = "eIF3_subunitE"
    eIF3_subunitF = "eIF3_subunitF"
    eIF3_subunitG = "eIF3_subunitG"
    eIF3_subunitH = "eIF3_subunitH"
    eIF3_subunitI = "eIF3_subunitI"
    eIF3_subunitJ = "eIF3_subunitJ"
    eIF3_subunitK = "eIF3_subunitK"
    eIF3_subunitL = "eIF3_subunitL"
    eIF3_subunitM = "eIF3_subunitM"

    eIF4F_4A = "eIF4F_4A"
    eIF4F_4G = "eIF4F_4G"
    eIF4F_4E = "eIF4F_4E"

    eIF4B = "eIF4B"
    eIF5B = "eIF5B"
    eIF5  = "eIF5"

    #!Bacterial
    IF1 = "IF1"
    IF2 = "IF2"
    IF3 = "IF3"

    #!Archaeal
    aIF_1A       = "aIF1A"
    aIF_2_alpha  = "aIF2_alpha"
    aIF_2_beta   = "aIF2_beta"
    aIF_2_gamma  = "aIF2_gamma"
    aIF_2B_alpha = "aIF2B_alpha"
    aIF_2B_beta  = "aIF2B_beta"
    aIF_2B_delta = "aIF2B_delta"
    aIF5A        = "aIF5A"
    aIF5B        = "aIF5B"

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


def get_polymer_color(polymer_class:str):
    polyix      = list(map(lambda x: x.value, (PolymerClass)))
    class_index = polyix.index(polymer_class) % len(polyix)

    if polymer_class == None:
        return "gray"
    return CHIMERAX_COLORS[(class_index % len(CHIMERAX_COLORS))][1]

def ribosome_representation(session, structure: AtomicStructure):
    for _ in list(PolymerClass):
        print(_.name)

    from chimerax.core.commands import run
    from chimerax.core.colors import hex_color
    from chimerax.atomic import Residue, Atom, Chain

    rcsb_id = str(structure.name).upper().split('.')[0] # <-- the structure gets opened with the basename ex "(5AFI.cif)" 

    run(session, "hide #2")

    with open(os.path.join(RIBETL_DATA, rcsb_id, "{}.json".format(rcsb_id)), "r") as f:
        profile = json.load(f)

    polymers = {}
    polymer_chains = [
        *profile["proteins"],
        *profile["rnas"],
        *profile["other_polymers"],
    ]

    [ polymers.update(x) for x in [{chain["auth_asym_id"]: chain} for chain in polymer_chains] ]

    c: Chain
    for c in structure.chains:

        aaid = c.chain_id
        polyclass = None

        if len(polymers[aaid]["nomenclature"]) < 1:
            continue
        else:
            polyclass = polymers[aaid]["nomenclature"][0]

        if polymers[aaid]["entity_poly_polymer_type"] == "RNA":
            run(session, "surf /{}".format(aaid))
            run(session, "color /{} gray".format(aaid))
        else:
            run(session, "show /{} cartoon".format(aaid, get_polymer_color(polyclass)))
            run(session, "color /{} {}".format(aaid, get_polymer_color(polyclass)))

    run(session, "set bgColor white")
    run(session, "graphics silhouettes true width 1")
    run(session, "light soft")

def produce_and_save_movie(session, target:str):
    print("GOT TARGET", target)
    RCSB_ID = target
    
    run(session, "open /home/rtviii/dev/RIBETL_DATA/{}/{}.cif".format(RCSB_ID, RCSB_ID))
    run(session, "sym #1 assembly 1") # take only one assembly if multiple are available
    run(session, "ribrep #2")
    run(session, "movie record")
    run(session, "turn y 2 180")
    run(session, "wait 180")
    run(session, "movie encode /home/rtviii/dev/riboxyz/chimerax/movies/{}.mp4".format(RCSB_ID))
    run(session, "close all")


def register_ribrepr_command(logger):
    from chimerax.core.commands import CmdDesc, register
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom

    desc = CmdDesc(
        required           = [("structure", AtomicStructureArg)],
        required_arguments = ["structure"],
        synopsis           = "representation ",
    )
    register("ribrep", desc, ribosome_representation, logger=logger)

def register_movie_command(logger):
    from chimerax.core.commands import CmdDesc, register, StringArg
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom
    desc = CmdDesc(
        required           = [("target", StringArg)],
        # required_arguments = ["structure"],
        synopsis           = "target ",
    )
    register("ribmovie", desc, produce_and_save_movie, logger=logger)


register_ribrepr_command(session.logger)
register_movie_command(session.logger)


