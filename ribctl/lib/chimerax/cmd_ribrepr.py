import enum
import json
import os
from chimerax.core.commands import register, CmdDesc
from chimerax.atomic import Structure, AtomicStructure, Chain
from chimerax.core.commands import run, runscript
from chimerax.core.commands import CmdDesc, register, StringArg
from chimerax.core.colors import hex_color

RIBETL_DATA = os.environ.get("RIBETL_DATA")




CytosolicProteinsColorScheme = {
    "uS2": "#60a5fa",
    "uS3": "#7dd3fc",
    "uS4": "#93c5fd",
    "uS5": "#38bdf8",
    "uS7": "#0ea5e9",
    "uS8": "#0284c7",
    "uS9": "#7dd3fc",
    "uS10": "#93c5fd",
    "uS11": "#38bdf8",
    "uS12": "#60a5fa",
    "uS13": "#3b82f6",
    "uS14": "#2563eb",
    "uS15": "#0ea5e9",
    "uS17": "#7dd3fc",
    "uS19": "#93c5fd",
    "bS1": "#6ee7b7",
    "bS6": "#34d399",
    "bS16": "#a7f3d0",
    "bS18": "#6ee7b7",
    "bS20": "#34d399",
    "bS21": "#10b981",
    "bTHX": "#059669",
    "eS1": "#c084fc",
    "eS4": "#d8b4fe",
    "eS6": "#e9d5ff",
    "eS7": "#a855f7",
    "eS8": "#9333ea",
    "eS10": "#c084fc",
    "eS12": "#d8b4fe",
    "eS17": "#a855f7",
    "eS19": "#c084fc",
    "eS21": "#d8b4fe",
    "eS24": "#e9d5ff",
    "eS25": "#a855f7",
    "eS26": "#9333ea",
    "eS27": "#c084fc",
    "eS28": "#d8b4fe",
    "eS30": "#e9d5ff",
    "eS31": "#a855f7",
    "RACK1": "#9333ea",
    "uL1": "#f59e0b",
    "uL2": "#dc2626",
    "uL3": "#fbbf24",
    "uL4": "#b45309",
    "uL5": "#e11d48",
    "uL6": "#d97706",
    "uL10": "#92400e",
    "uL11": "#fef9c3",
    "uL13": "#9a3412",
    "uL14": "#c2410c",
    "uL15": "#fb923c",
    "uL16": "#be123c",
    "uL18": "#fde047",
    "uL22": "#92400e",
    "uL23": "#ff8fab",
    "uL24": "#b45309",
    "uL29": "#f97316",
    "uL30": "#ea580c",
    "bL9": "#854d0e",
    "bL12": "#facc15",
    "bL17": "#d97706",
    "bL19": "#991b1b",
    "bL20": "#fb923c",
    "bL21": "#ef4444",
    "bL25": "#78350f",
    "bL27": "#fb7185",
    "bL28": "#92400e",
    "bL31": "#713f12",
    "bL32": "#fef08a",
    "bL33": "#ea580c",
    "bL34": "#b91c1c",
    "bL35": "#9a3412",
    "bL36": "#f59e0b",
    "eL6": "#fef3c7",
    "eL8": "#b45309",
    "eL13": "#92400e",
    "eL14": "#dc2626",
    "eL15": "#fbbf24",
    "eL18": "#9f1239",
    "eL19": "#f97316",
    "eL20": "#fdba74",
    "eL21": "#92400e",
    "eL22": "#ff8fab",
    "eL24": "#c2410c",
    "eL27": "#fde047",
    "eL28": "#e11d48",
    "eL29": "#fef9c3",
    "eL30": "#92400e",
    "eL31": "#fb923c",
    "eL32": "#b91c1c",
    "eL33": "#ea580c",
    "eL34": "#f59e0b",
    "eL36": "#854d0e",
    "eL37": "#fef08a",
    "eL38": "#dc2626",
    "eL39": "#9a3412",
    "eL40": "#fbbf24",
    "eL41": "#9f1239",
    "eL42": "#fb923c",
    "eL43": "#78350f",
    "P1P2": "#e11d48",
}
MitochondrialProteinColorScheme = {
    "uS2m":"#67e8f9",
    "uS3m":"#22d3ee",
    "uS4m":"#7dd3fc",
    "uS5m":"#38bdf8",
    "uS7m":"#0ea5e9",
    "uS8m":"#0284c7",
    "uS9m":"#2dd4bf",
    "uS10m":"#14b8a6",
    "uS11m":"#0d9488",
    "uS12m":"#5eead4",
    "uS13m":"#a5f3fc",
    "uS14m":"#7dd3fc",
    "uS15m":"#38bdf8",
    "uS17m":"#0ea5e9",
    "uS19m":"#0284c7",
    "bS1m":"#a3e635",
    "bS6m":"#84cc16",
    "bS16m":"#bef264",
    "bS18m":"#d9f99d",
    "bS21m":"#65a30d",
    "mS22":"#fcaec8",
    "mS23":"#fda4af",
    "mS25":"#fba5b5",
    "mS26":"#fecaca",
    "mS27":"#fecdd3",
    "mS29":"#ffd4d4",
    "mS31":"#ffe4e6",
    "mS33":"#fce7f3",
    "mS34":"#fdf2f8",
    "mS35":"#fbb4c7",
    "mS37":"#fda4af",
    "mS38":"#fecdd3",
    "mS39":"#fee2e2",
    "mS40":"#fecaca",
    "mS41":"#fee2e2",
    "mS42":"#fef2f2",
    "mS43":"#fbb4c7",
    "mS44":"#fda4af",
    "mS45":"#fecdd3",
    "mS46":"#ffb3c1",
    "mS47":"#ff8fab",
    "uL1m":"#fda4af",
    "uL2m":"#fecdd3",
    "uL3m":"#fee2e2",
    "uL4m":"#fecaca",
    "uL5m":"#fee2e2",
    "uL6m":"#fef2f2",
    "uL10m":"#ff8fab",
    "uL11m":"#ffb3c1",
    "uL13m":"#fbb4c7",
    "uL14m":"#fda4af",
    "uL15m":"#fecdd3",
    "uL16m":"#fee2e2",
    "uL18m":"#fecaca",
    "uL22m":"#fee2e2",
    "uL23m":"#fef2f2",
    "uL24m":"#ffb3c1",
    "uL29m":"#fba4af",
    "uL30m":"#fecdd3",
    "bL9m":"#fef9c3",
    "bL12m":"#fef08a",
    "bL17m":"#fde047",
    "bL19m":"#facc15",
    "bL20m":"#eab308",
    "bL21m":"#fbbf24",
    "bL27m":"#fcd34d",
    "bL28m":"#fde68a",
    "bL31m":"#fef3c7",
    "bL32m":"#fffbeb",
    "bL33m":"#fde047",
    "bL34m":"#facc15",
    "bL35m":"#eab308",
    "bL36m":"#fbbf24",
    "mL37":"#fae8ff",
    "mL38":"#f3e8ff",
    "mL39":"#f0abfc",
    "mL40":"#e879f9",
    "mL41":"#d946ef",
    "mL42":"#c084fc",
    "mL43":"#a855f7",
    "mL44":"#9333ea",
    "mL45":"#fae8ff",
    "mL46":"#f3e8ff",
    "mL48":"#f0abfc",
    "mL49":"#e879f9",
    "mL50":"#d946ef",
    "mL51":"#c084fc",
    "mL52":"#a855f7",
    "mL53":"#9333ea",
    "mL54":"#fae8ff",
    "mL57":"#f3e8ff",
    "mL58":"#f0abfc",
    "mL59":"#e879f9",
    "mL60":"#d946ef",
    "mL61":"#c084fc",
    "mL62":"#a855f7",
    "mL63":"#9333ea",
    "mL64":"#fae8ff",
    "mL65":"#f3e8ff",
    "mL66":"#f0abfc",
    "mL67":"#e879f9",
}
RNAColorScheme = {
    "5SrRNA"   :"#e2e8f0",
    "16SrRNA"  :"#bfdbfe",
    "23SrRNA"  :"#93c5fd",
    "25SrRNA"  :"#f8fafc",
    "5.8SrRNA" :"#94a3b8",
    "18SrRNA"  :"#d1d5db",
    "28SrRNA"  :"#fafafa",
    "mt12SrRNA":"#e5e7eb",
    "mt16SrRNA":"#f1f5f9",
    "tRNA"     :"#c084fc",
}
FactorsColorScheme = {
    "eEF1A":"#fed7aa",
    "eEF1B":"#fdba74",
    "eFSec":"#fb923c",
    "eEF2" :"#fcd34d",
    "mtEF4":"#fde047",
    "eIF5A":"#fef9c3",
    "eEF3" :"#fef3c7",
    "EF-Tu":"#fecaca",
    "EF-Ts":"#fecdd3",
    "SelB" :"#fee2e2",
    "EF-G" :"#fef2f2",
    "EF4"  :"#ffe4e6",
    "EF-P" :"#fce7f3",
    "Tet_O":"#fed7aa",
    "Tet_M":"#fdba74",
    "RelA" :"#fb923c",
    "BipA" :"#fcd34d",
    "aEF1A":"#fef3c7",
    "aEF2" :"#fde68a",
}
InitiationFactorsColorScheme = {
    "eIF1"         : "#fef9c3",
    "eIF1A"        : "#fef08a",
    "eIF2_alpha"   : "#fde047",
    "eIF2_beta"    : "#facc15",
    "eIF2_gamma"   : "#eab308",
    "eIF2B_alpha"  : "#fef9c3",
    "eIF2B_beta"   : "#fef08a",
    "eIF2B_gamma"  : "#fde047",
    "eIF2B_delta"  : "#facc15",
    "eIF2B_epsilon": "#eab308",
    "eIF3_subunitA": "#fef3c7",
    "eIF3_subunitB": "#fde68a",
    "eIF3_subunitC": "#fcd34d",
    "eIF3_subunitD": "#fbbf24",
    "eIF3_subunitE": "#f59e0b",
    "eIF3_subunitF": "#fef9c3",
    "eIF3_subunitG": "#fef08a",
    "eIF3_subunitH": "#fde047",
    "eIF3_subunitI": "#facc15",
    "eIF3_subunitJ": "#eab308",
    "eIF3_subunitK": "#fef3c7",
    "eIF3_subunitL": "#fde68a",
    "eIF3_subunitM": "#fcd34d",
    "eIF4F_4A"     : "#fef9c3",
    "eIF4F_4G"     : "#fef08a",
    "eIF4F_4E"     : "#fde047",
    "eIF4B"        : "#facc15",
    "eIF5B"        : "#eab308",
    "eIF5"         : "#fef9c3",
    "IF1"          : "#f8fafc",
    "IF2"          : "#f1f5f9",
    "IF3"          : "#e2e8f0",
    "aIF1A"        : "#fef9c3",
    "aIF2_alpha"   : "#fef08a",
    "aIF2_beta"    : "#fde047",
    "aIF2_gamma"   : "#facc15",
    "aIF2B_alpha"  : "#eab308",
    "aIF2B_beta"   : "#fef9c3",
    "aIF2B_delta"  : "#fef08a",
    "aIF5A"        : "#fde047",
    "aIF5B"        : "#facc15",
}



POLYMER_COLORS = {
    **CytosolicProteinsColorScheme,
    **MitochondrialProteinColorScheme,
    **RNAColorScheme,
    **FactorsColorScheme,
    **InitiationFactorsColorScheme
}












def get_polymer_color(polymer_class: str) -> str:
    """Get the hex color for a polymer class, defaulting to gray if not found."""
    return POLYMER_COLORS.get(polymer_class, "#808080")

def ribosome_representation(session, structure: AtomicStructure):
    from chimerax.core.commands import run
    from chimerax.core.colors import hex_color
    from chimerax.atomic import Residue, Atom, Chain

    rcsb_id = str(structure.name).upper().split('.')[0]
    run(session, "set bgColor white")
    run(session, "sym #1 assembly 1")
    run(session, "hide #2")
    run(session, "hide") 

    profile_path = os.path.join(RIBETL_DATA, rcsb_id, f"{rcsb_id}.json")
    with open(profile_path, "r") as f:
        profile = json.load(f)

    # Create chain-to-polymer mapping
    polymers = {}
    polymer_chains = [*profile["proteins"], *profile["rnas"], *profile["other_polymers"]]
    for chain in polymer_chains:
        polymers[chain["auth_asym_id"]] = chain

    # Process each chain
    for c in structure.chains:
        chain_id = c.chain_id
        if chain_id not in polymers:
            continue

        chain_info = polymers[chain_id]
        polyclass = chain_info["nomenclature"][0] if chain_info["nomenclature"] else None
        
        # Get color based on polymer class
        color = get_polymer_color(polyclass)
        
        if chain_info["entity_poly_polymer_type"] == "RNA":
            run(session, f"surf /{chain_id}")
            run(session, f"transp /{chain_id} 100")
            run(session, f"color /{chain_id} {color}")
        else:
            # For proteins, use a single show command with cartoon style
            run(session, f"show /{chain_id} cartoon")
            run(session, f"color /{chain_id} {color}")

    # Final styling
    run(session, "graphics silhouettes true width 1")
    run(session, "light soft")

def register_ribrepr_command(logger):

    from chimerax.core.commands import CmdDesc, register
    from chimerax.atomic import AtomicStructureArg, Chain, Residue, Atom
    desc = CmdDesc(
        required=[("structure", AtomicStructureArg)],
        required_arguments=["structure"],
        synopsis="Apply ribosome representation with custom coloring"
    )
    register("ribrep", desc, ribosome_representation, logger=logger)

register_ribrepr_command(session.logger)