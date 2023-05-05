import argparse
import json
import os
import sys
from pymol import cmd

RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))

def render_thumbnail(pdbid: str):

    pdbid = pdbid.upper()
    thumbnail_path = os.path.join(RIBETL_DATA, pdbid, f"_ray_{pdbid}.png")

    if os.path.exists(thumbnail_path):
        print(f"Thumbnail already exists: {thumbnail_path}")
        sys.exit(0)

    else:
        cmd.load(os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif"))
        cmd.reset()
        cmd.spectrum('chain')
        cmd.ray(500, 500)
        cmd.png(thumbnail_path)
        print('Saved {}'.format(thumbnail_path))


if __name__ == "__main__":
    from pymol import cmd, util

    parser = argparse.ArgumentParser(
        description='Generate a ribosome thumbnail image.')
    parser.add_argument("-s", "--structure", type=str,
                        required=True, help="RCSB ID of structure to process")
    args = parser.parse_args()
    pdbid = args.structure.upper()
    thumbnail_path = os.path.join(RIBETL_DATA, pdbid, f"_ray_{pdbid}.png")

    if os.path.exists(thumbnail_path):
        print(f"Thumbnail already exists: {thumbnail_path}")
        sys.exit(0)

    else:
        cmd.load(os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif"))
        cmd.reset()
        cmd.spectrum('chain')
        cmd.ray(500, 500)
        cmd.png(thumbnail_path)
        print('Saved {}'.format(thumbnail_path))

colormap__LSU_Proteins = {
                          "uL1": "lawrencium",
                          "uL2": "lead",
                          "uL3": "lithium",
                          "uL4": "lutetium",
                          "uL5": "magnesium",
                          "uL6": "manganese",
                          "eL6": "meitnerium",
                          "eL8": "mendelevium",
                          "bL9": "mercury",
                          "uL10": "molybdenum",
                          "uL11": "neodymium",
                          "bL12": "neon",
                          "uL13": "neptunium",
                          "eL13": "nickel",
                          "uL14": "niobium",
                          "eL14": "nitrogen",
                          "uL15": "nobelium",
                          "eL15": "osmium",
                          "uL16": "oxygen",
                          "bL17": "palladium",
                          "uL18": "phosphorus",
                          "eL18": "platinum",
                          "bL19": "plutonium",
                          "eL19": "polonium",
                          "bL20": "potassium",
                          "eL20": "praseodymium",
                          "bL21": "promethium",
                          "eL21": "protactinium",
                          "uL22": "radium",
                          "eL22": "radon",
                          "uL23": "rhenium",
                          "uL24": "rhodium",
                          "eL24": "rubidium",
                          "bL25": "ruthenium",
                          "bL27": "rutherfordium",
                          "eL27": "samarium",
                          "bL28": "scandium",
                          "eL28": "seaborgium",
                          "uL29": "selenium",
                          "eL29": "silicon",
                          "uL30": "silver",
                          "eL30": "sodium",
                          "bL31": "strontium",
                          "eL31": "sulfur",
                          "bL32": "tantalum",
                          "eL32": "technetium",
                          "bL33": "tellurium",
                          "eL33": "terbium",
                          "bL34": "thallium",
                          "eL34": "thorium",
                          "bL35": "thulium",
                          "bL36": "tin",
                          "eL36": "titanium",
                          "eL37": "tungsten",
                          "eL38": "uranium",
                          "eL39": "vanadium",
                          "eL40": "xenon",
                          "eL41": "ytterbium",
                          "eL42": "yttrium",
                          "eL43": "zinc",
                          "P1/P2": "zirconium"
                          
                          }

colormap__SSU_Proteins = {
                     "bS1": "actinium",
                     "eS1": "aluminum",
                     "uS2": "americium",
                     "uS3": "antimony",
                     "uS4": "argon",
                     "eS4": "arsenic",
                     "uS5": "astatine",
                     "bS6": "barium",
                     "eS6": "berkelium",
                     "uS7": "beryllium",
                     "eS7": "bismuth",
                     "uS8": "bohrium",
                     "eS8": "boron",
                     "uS9": "bromine",
                     "uS10": "cadmium",
                     "eS10": "calcium",
                     "uS11": "californium",
                     "uS12": "carbon",
                     "eS12": "cerium",
                     "uS13": "cesium",
                     "uS14": "chlorine",
                     "uS15": "chromium",
                     "bS16": "cobalt",
                     "uS17": "copper",
                     "eS17": "curium",
                     "bS18": "deuterium",
                     "uS19": "dubnium",
                     "eS19": "dysprosium",
                     "bS20": "einsteinium",
                     "bS21": "erbium",
                     "bTHX": "europium",
                     "eS21": "fermium",
                     "eS24": "fluorine",
                     "eS25": "francium",
                     "eS26": "gadolinium",
                     "eS27": "gallium",
                     "eS28": "germanium",
                     "eS30": "gold",
                     "eS31": "hafnium",
                     }