
import argparse
import json
import os
import sys
from pymol import cmd
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))


def sload(pdbid: str):

    pdbid = pdbid.upper()
    RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
    path = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.cif")
    cmd.load(path)


def render_thumbnail(pdbid: str):

    from pymol import cmd

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

def by_rna(pdbid: str):
    pdbid = pdbid.upper()
    RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
    path = os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.json")

    with open(path, 'rb') as infile:
        profile = json.load(infile)

    if profile['rnas'] != None:
        for rna in profile['rnas']:
            # cmd.color('white', f"chain {rna['auth_asym_id']}")
            cmd.hide('everything', f"chain {rna['auth_asym_id']}")
            cmd.show('surface', f"chain {rna['auth_asym_id']}")
            cmd.color('white', f"chain {rna['auth_asym_id']}")
            cmd.set('transparency', 0.1, f"chain {rna['auth_asym_id']}")

    for protein in profile['proteins']:

        if len(protein['nomenclature']) != 0:
            if protein['nomenclature'][0] in colormap__LSU_Proteins:
                CLR = colormap__LSU_Proteins[protein['nomenclature'][0]]
            elif protein['nomenclature'][0] in colormap__SSU_Proteins:
                CLR = colormap__SSU_Proteins[protein['nomenclature'][0]]
        else:
            CLR = 'blue'

        prot_tmp = f"protein_tmp_{ protein['auth_asym_id'] }"

        cmd.create(prot_tmp, f"chain {protein[ 'auth_asym_id' ]}")
        cmd.remove(f"chain {protein[ 'auth_asym_id' ]} and m. {pdbid}")
        cmd.hide ('everything', prot_tmp)
        cmd.show ('surface', prot_tmp)
        cmd.show ('sticks', prot_tmp)
        cmd.show ('cartoon', prot_tmp)
        cmd.show ('lines', prot_tmp)
        cmd.set('transparency', 0.75, prot_tmp)
        cmd.color(CLR , prot_tmp)


cmd.extend("sload", sload)
cmd.extend("by_rna", by_rna)

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