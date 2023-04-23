from functools import reduce
import os
from pprint import pprint
import sys
from fuzzysearch import find_near_matches
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom
from api.ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets, obtain_assets, obtain_assets_threadpool
from api.ribctl.lib.utils import open_structure

RIBETL_DATA = os.environ.get('RIBETL_DATA')
# ※ ---------------------------- 23/25/28SrRNA PTC residue locations ---------------------------- ※
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "subseq_6": "AAGACCC",
    "subseq_8": "GGAUAAC",
    "subseq_9": "GAGCUGGGUUUA",
    "e": {
        "site_6": [2400, 2401, 2402, 2403, 2404, 2405, 2406, 2407],
        "site_8": [2811, 2812, 2813, 2814, 2815, 2816, 2817, 2818],
        "site_9": [2941, 2942, 2943, 2944, 2945, 2946, 2947, 2948, 2949, 2950, 2951, 2952, 2953]
    },
    'b': {
        "site_6": [2059, 2060, 2061,
                   2062, 2063, 2064,
                   2065],
        "site_8": [2446, 2447, 2448,
                   2449, 2450, 2451, 2452],
        "site_9": [2576, 2577, 2578, 2579,
                   2580, 2581, 2582, 2583,
                   2584, 2585, 2586, 2587],
    }
}



# RCSB_ID = '3J92'
RCSB_ID = sys.argv[1].upper()


def pick_match(ms, rna_length: int):
    if len(ms) == 0:
        print("No matches found!!")
        exit(1)

    """Pick the match that is closest to the 3' end of the rRNA."""
    best = None
    farthest_dist = 0

    if len(ms) > 1:
        for m in ms:
            if rna_length - ((m.start + m.end) / 2) < farthest_dist:
                best = m
        return best
    else:
        return ms[0]

# def ptc_by_alignment(rcsb_id: str, organism:int):
def ptc_by_alignment():
    R = RibosomeAssets(RCSB_ID)
    profile = R.profile()
    
    
    struct_profile: Structure = R.biopython_structure()
    [ strand_target,chain_id, rna_type] = R.get_LSU_rRNA()
    

    
    profile.assembly_map[0]
    
    if chain_id in struct_profile.child_dict[0].child_dict:
        chain3d: Chain = struct_profile.child_dict[0].child_dict[chain_id]
    else:
        chain3d: Chain = struct_profile.child_dict[1].child_dict[chain_id]

    ress:list[Residue] = chain3d.child_list
    ress_sanitized: list[Residue] = [*filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"


    raw_seq         = reduce(lambda x,y: x + y.resname, ress_sanitized,'')
    matches         = find_near_matches(DORIS_ET_AL["subseq_9"], raw_seq, max_l_dist=1)
    m0              = pick_match(matches, len(raw_seq))
    sought_residues = [ ress_sanitized[i] for i in list(range(m0.start, m0.end))]

    PTC_MARKERKS_RAW = {}
    for res in sought_residues:
        if res.id[1] not in PTC_MARKERKS_RAW:
            PTC_MARKERKS_RAW[res.id[1]] = {}
        atom:Atom
        for atom in res.child_list:
            atom_name = atom.name
            atom_coords = atom.get_coord()
            PTC_MARKERKS_RAW[res.id[1]][atom_name] = list(map(lambda x: float(x), list(atom_coords)))

    PTC_MARKERKS_RAW = {
        chain_id: PTC_MARKERKS_RAW
    }

    markers_dir = os.path.join(RIBETL_DATA, "PTC_MARKERS_RAW")
    outfile     = os.path.join(markers_dir, f"{rcsb_id.upper()}_PTC_MARKERS_RAW.json")

    with open(outfile, 'w') as outf:
        json.dump(PTC_MARKERKS_RAW, outf, indent=4)
        print("Saved {} successfully.".format(outfile))


R       = RibosomeAssets(RCSB_ID)
profile = R.profile()
rna = R.get_LSU_rRNA()
seq, auth_asym_id, rna_type = rna.entity_poly_seq_one_letter_code_can, rna.auth_asym_id, rna.nomenclature[0] if len(rna.nomenclature) > 0 else None


# pprint(profile.assembly_map[0].dict())
# pprint(profile.assembly_map[1].dict())
# pprint(auth_asym_id)


print(R.get_rna_by_nomclass("23SrRNA"))

# def extract_ptc_coordinates():
#     ...



# ./6HCM/polymer_1_nascent_chain.json
# ./7A5I/polymer_Y2_nascent_chain.json
# ./6HCF/polymer_1_nascent_chain.json
# ./3JAJ/polymer_2_nascent_chain.json
# ./4V5H/polymer_AZ_nascent_chain.json
# ./3J7Z/polymer_a_nascent_chain.json
# ./5LZW/polymer_1_nascent_chain.json
# ./7A5G/polymer_Y2_nascent_chain.json
# ./5LZZ/polymer_1_nascent_chain.json
# ./6SGC/polymer_XX_nascent_chain.json
# ./7A5F/polymer_Y2_nascent_chain.json
# ./6OLI/polymer_y_nascent_chain.json
# ./5LZX/polymer_1_nascent_chain.json
# ./6HCJ/polymer_1_nascent_chain.json
# ./5LZT/polymer_1_nascent_chain.json
# ./6HCQ/polymer_1_nascent_chain.json
# ./7A5K/polymer_Y2_nascent_chain.json
# ./3J92/polymer_1_nascent_chain.json
# ./5NWY/polymer_s_nascent_chain.json
# ./7A5H/polymer_Y2_nascent_chain.json
# ./6XA1/polymer_NC_nascent_chain.json
# ./4V6M/polymer_AV_nascent_chain.json
# ./6W6L/polymer_y_nascent_chain.json
# ./5LZV/polymer_1_nascent_chain.json
# ./7A5J/polymer_Y2_nascent_chain.json