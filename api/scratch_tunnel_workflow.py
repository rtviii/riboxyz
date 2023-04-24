from functools import reduce
import os
from pprint import pprint
import sys
from fuzzysearch import find_near_matches
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Atom import Atom
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_binding_site import AMINO_ACIDS, NUCLEOTIDES, ResidueSummary
from api.ribctl.lib.types.types_ribosome import RibosomeStructure

RIBETL_DATA = os.environ.get('RIBETL_DATA')

# ※ ---------------------------- 23/25/28SrRNA PTC residue locations ---------------------------- ※
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "subseq_6": "AAGACCC",
    "subseq_8": "GGAUAAC",
    "subseq_9": "GAGCUGGGUUUA",
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
# ※ --------------------------------------------------------------------------------------------- ※

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


def ptc_residues_via_alignment(rcsb_id: str, assembly_id: int = 0) -> tuple[list[Residue], str]:
    R = RibosomeAssets(rcsb_id)
    # R.profile()

    struct_profile: Structure = R.biopython_structure()
    rna = R.get_LSU_rRNA(assembly_id)
    auth_asym_id = rna.auth_asym_id

    chain3d: Chain = struct_profile.child_dict[assembly_id].child_dict[auth_asym_id]
    ress: list[Residue] = chain3d.child_list
    ress_sanitized: list[Residue] = [
        *filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"

    # create sequence string from residues
    raw_seq = reduce(lambda x, y: x + y.resname, ress_sanitized, '')
    # find a subsequences that matches doris with max levenstein distance of 1

    matches = find_near_matches(DORIS_ET_AL["subseq_9"], raw_seq, max_l_dist=1)
    m0 = pick_match(matches, len(raw_seq))
    PTC_residues = [ress_sanitized[i] for i in list(range(m0.start, m0.end))]

    return PTC_residues, auth_asym_id


def ptc_coordinates(
    reslist: list[Residue],
    auth_asym_id: str
) -> dict:
    """
    chain:
        residue_id:
            atom_name: [x, y, z]
    """

    ptc_coordinates = {}
    for res in reslist:
        if res.id[1] not in ptc_coordinates:
            ptc_coordinates[res.id[1]] = {}

        atom: Atom
        for atom in res.child_list:
            atom_name = atom.name
            atom_coords = atom.get_coord()
            ptc_coordinates[res.id[1]][atom_name] = list(
                map(lambda x: float(x), list(atom_coords)))

    ptc_coordinates = {
        auth_asym_id: ptc_coordinates
    }

    return ptc_coordinates


def ptc_midpoint(reslist: list[Residue], auth_asym_id: str) -> tuple[float, float, float]:
    """
    chain:
        residue_id:
            atom_name: [x, y, z]

    """

    ptc_coord_dict = ptc_coordinates(reslist, auth_asym_id)

    lsu_rna = [*ptc_coord_dict.keys()][0]
    if lsu_rna == None:
        raise Exception("Could not identify chain")

    subseq_residues = ptc_coord_dict[lsu_rna]
    subseq_length = len(subseq_residues)

    if "O4'" in [*subseq_residues.values()][subseq_length - 2]:
        # pre  last residue of the comb
        U_end_pos = [*subseq_residues.values()][subseq_length - 2]["O4'"]
    else:
        # pre  last residue of the comb
        U_end_pos = [*subseq_residues.values()][subseq_length - 2]["C4"]

    if "O4'" in [*subseq_residues.values()][0]:
        # first residue of the comb
        U_start_pos = [*subseq_residues.values()][0]["O4'"]
    else:
        # first residue of the comb
        U_start_pos = [*subseq_residues.values()][0]["C4"]

    midpoint = (
        (U_end_pos[0] + U_start_pos[0]) / 2,
        (U_end_pos[1] + U_start_pos[1]) / 2,
        (U_end_pos[2] + U_start_pos[2]) / 2,
    )

    return midpoint

def residue_is_canonical(res: Residue | ResidueSummary) -> bool:
    return res.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]


def tunnel_obstructions(rcsb_id:str, ptc_midpoint:tuple[float,float,float], radius:int=30)->list[dict]:
    
    R = RibosomeAssets(RCSB_ID)
    structure, profile = R.get_struct_and_profile()

    neigbor_search = NeighborSearch(list(structure.get_atoms()))
    nbr_residues   = []

    nbr_chains = neigbor_search.search(ptc_midpoint, radius, level='C')
    print("nbr chains")
    for c in nbr_chains:
        print(c.id)
    
    nbr_residues.extend(neigbor_search.search(ptc_midpoint, radius, level='R'))
    nbr_residues     = list(set([* map(ResidueSummary.from_biopython_residue, nbr_residues)]))
    foreign_residues = []


    for N in nbr_residues:
         parent = R.get_chain_by_auth_asym_id(N.parent_auth_asym_id)
         if parent is not None:
             print(parent.auth_asym_id,type(parent), residue_is_canonical(N))
         else:
             raise LookupError("Parent chain must exist in structure. Something went wrong (with the data, probably)")
    print("Inspected midpoint ", ptc_midpoint, " with radius", radius)
         



    # print("Looking for stuff aroudn midpoint ", midpoint)
    # foreign_residues = list(filter(lambda res: res.resname not in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues))
    # for F in foreign_residues:
    #     print(F)

    

    return nbr_residues


idlist = [
    "6HCM",
    "7A5I",
    "6HCF",
    "3JAJ",
    "4V5H",
    "3J7Z",
    "5LZW",
    "7A5G",
    "5LZZ",
    "6SGC",
    "7A5F",
    "6OLI",
    "5LZX",
    "6HCJ",
    "5LZT",
    "6HCQ",
    "7A5K",
    "3J92",
    "5NWY",
    "7A5H",
    "6XA1",
    "4V6M",
    "6W6L",
    "5LZV",
    "7A5J"
]

ptcres, auth_asym_id = ptc_residues_via_alignment(RCSB_ID, 0)
midpoint             = ptc_midpoint(ptcres, auth_asym_id)
tunnel_obstructions(RCSB_ID, midpoint)


# for rcsb in idlist:
#     try:
#         ptcres, auth_asym_id = ptc_residues_via_alignment(rcsb, 0)
#         print(rcsb, ptc_midpoint(ptcres, auth_asym_id))
#     except:
#         ...
