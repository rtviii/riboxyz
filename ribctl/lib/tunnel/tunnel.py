import math
from pprint import pprint


from ribctl.lib.msalib import (
    util__backwards_match,
    util__forwards_match,
)
from ribctl.lib.types.types_ribosome import PolymericFactor, Protein, ProteinClass
from ribctl.lib.types.types_binding_site import AMINO_ACIDS, NUCLEOTIDES, ResidueSummary
from ribctl.lib.mod_transpose_bsites import SeqMatch

from Bio.PDB.Atom import Atom
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from functools import reduce
import numpy as np
import os

from scripts.prd_to_biopython import msa_dict, msaclass_extend_temp



# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "SITE_6": "AAGACCC",
    "SITE_8": "GGAUAAC",
    "SITE_9": "GAGCUGGGUUUA",
}


def pick_match(matches, rna_length: int):
    if len(matches) == 0:
        return None
    """Pick the match that is closest to the 3' end of the rRNA."""
    best = None
    farthest_dist = 0

    if len(matches) > 1:
        for m in matches:
            if rna_length - ((m.start + m.end) / 2) < farthest_dist:
                best = m
        return best
    else:
        return matches[0]


def residue_labels(res: Residue | ResidueSummary) -> bool:
    return res.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]


def ptc_fuzzyfind_subseq_in_chain(biopython_struct, auth_asym_id:str, assembly_id:int=0):
    chain3d: Chain = biopython_struct.child_dict[assembly_id].child_dict[auth_asym_id]
    ress: list[Residue] = chain3d.child_list
    ress_sanitized: list[Residue] = [
        *filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "-", "PSU"], ress)
    ]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"

    # create sequence string from residues
    raw_seq = reduce(lambda x, y: x + y.resname, ress_sanitized, "")
    # find a subsequences that matches doris with max levenstein distance of 1

    # print("Raw seq is ", raw_seq)
    from fuzzysearch import find_near_matches

    match9         = pick_match(find_near_matches(DORIS_ET_AL["SITE_9"], raw_seq, max_l_dist=1), len(raw_seq))
    match8         = pick_match(find_near_matches(DORIS_ET_AL["SITE_8"], raw_seq, max_l_dist=1), len(raw_seq))
    match6         = pick_match(find_near_matches(DORIS_ET_AL["SITE_6"], raw_seq, max_l_dist=1), len(raw_seq))

    PTC_residues_9 = [ress_sanitized[i] for i in list(range(match9.start, match9.end))] if match9 else []
    PTC_residues_8 = [ress_sanitized[i] for i in list(range(match8.start, match8.end))] if match8 else []
    PTC_residues_6 = [ress_sanitized[i] for i in list(range(match6.start, match6.end))] if match6 else []

    print("Found {}, {}, {} in  ".format(match6,match8,match9), auth_asym_id, biopython_struct)
    return PTC_residues_6, PTC_residues_8, PTC_residues_9, auth_asym_id

def ptc_resdiues_get(rcsb_id: str, assembly_id: int = 0) -> tuple[list[Residue], str]:
    """
    @ rcsb_id - str
    @ assembly_id - int
    Given a bacterial structure (assertion unenforced):
     - obtain the LSU rRNA chain
     - sanitze the sequence (canonicalize non-canonical residues)
     - locate the conserved subsequence (due to Doris et al.)
     Returns the residue list for the ptc and the `auth_asym_id` of the rRNA chain
    """

    R = RibosomeAssets(rcsb_id)
    struct_profile = R.biopython_structure()
    try:
        rna = R.get_LSU_rRNA(assembly_id)
    except Exception as e:
        print("Unfortunately, ", e)
        rnas = R.profile().rnas
        potential_rnas = {}
        for _potential_lsurna in rnas:
            if len(_potential_lsurna.entity_poly_seq_one_letter_code_can) > 2500:
                potential_rnas = {
                    _potential_lsurna.rcsb_pdbx_description: _potential_lsurna,
                    **potential_rnas,
                }
        if len(potential_rnas.keys()) ==0:
            raise Exception("No 28S-like eukaryotic rRNAs in this struct")
        else:
            matches={}
            for p_rna in potential_rnas:
                print("processing ", p_rna)
                try:
                    m6,m8,m9, aaid = ptc_fuzzyfind_subseq_in_chain(struct_profile, potential_rnas[p_rna].auth_asym_id)
                    matches = {
                        **matches,
                        aaid: [m6,m8,m9]
                        }
                    if m9 != None:
                        rna = potential_rnas[p_rna]
                except Exception as e :
                    print(e)
                    ...



    auth_asym_id = rna.auth_asym_id

    chain3d: Chain = struct_profile.child_dict[assembly_id].child_dict[auth_asym_id]
    ress: list[Residue] = chain3d.child_list
    ress_sanitized: list[Residue] = [
        *filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)
    ]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"

    # create sequence string from residues
    raw_seq = reduce(lambda x, y: x + y.resname, ress_sanitized, "")
    # find a subsequences that matches doris with max levenstein distance of 1

    from fuzzysearch import find_near_matches

    matches = find_near_matches(DORIS_ET_AL["SITE_9"], raw_seq, max_l_dist=1)
    m0 = pick_match(matches, len(raw_seq))
    PTC_residues = [ress_sanitized[i] for i in list(range(m0.start, m0.end))]

    return PTC_residues, auth_asym_id


def ptc_residues_to_atom_coordinates(
    reslist: list[Residue], auth_asym_id: str
) -> dict[str, dict[str, list[float]]]:
    """
    Given a list of residues (representing the PTC), return a dictionary of the form:

    auth_asym_id:
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
                map(lambda x: float(x), list(atom_coords))
            )

    ptc_coordinates = {auth_asym_id: ptc_coordinates}

    return ptc_coordinates


def ptc_residues_calculate_midpoint(
    reslist: list[Residue], auth_asym_id: str
) -> list[float]:
    """
    auth_asym_id:
        residue_id:
            atom_name: [x, y, z]

    """

    ptc_coord_dict = ptc_residues_to_atom_coordinates(reslist, auth_asym_id)

    lsu_rna = [*ptc_coord_dict.keys()][0]
    if lsu_rna == None:
        raise Exception("Could not identify chain")

    subseq_residues = ptc_coord_dict[lsu_rna]
    subseq_length = len(subseq_residues)

    if "O4'" in [*subseq_residues.values()][subseq_length - 2]:
        # pre last residue of the comb
        U_end_pos = [*subseq_residues.values()][subseq_length - 2]["O4'"]
    else:
        # pre last residue of the comb
        U_end_pos = [*subseq_residues.values()][subseq_length - 2]["C4"]

    if "O4'" in [*subseq_residues.values()][0]:
        # first residue of the comb
        U_start_pos = [*subseq_residues.values()][0]["O4'"]
    else:
        # first residue of the comb
        U_start_pos = [*subseq_residues.values()][0]["C4"]

    midpoint = [
        (U_end_pos[0] + U_start_pos[0]) / 2,
        (U_end_pos[1] + U_start_pos[1]) / 2,
        (U_end_pos[2] + U_start_pos[2]) / 2,
    ]

    return midpoint


def tunnel_obstructions(
    rcsb_id: str, ptc_midpoint: tuple[float, float, float], radius: int = 30
) -> tuple[list[PolymericFactor], list[ResidueSummary]]:
    R = RibosomeAssets(rcsb_id)
    structure, profile = R.get_struct_and_profile()

    neigbor_search = NeighborSearch(list(structure.get_atoms()))
    nbr_residues = []

    nbr_residues.extend(neigbor_search.search(ptc_midpoint, radius, level="R"))
    nbr_residues = list(
        set([*map(ResidueSummary.from_biopython_residue, nbr_residues)])
    )

    foreign_nonpolys = []
    foregin_polymers = []

    for N in nbr_residues:
        if not residue_labels(N):
            foreign_nonpolys.append(N)
        parent_chain, chain_type = R.get_chain_by_auth_asym_id(N.parent_auth_asym_id)
        if (parent_chain, chain_type) != (None, None):
            if chain_type == "PolymericFactor":
                foregin_polymers.append(parent_chain)
        else:
            raise LookupError(
                "Parent chain must exist in structure. Something went wrong (with the data, probably)"
            )

    print("Inspected midpoint ", ptc_midpoint, " with radius", radius)
    return list(set(foregin_polymers)), list(set(foreign_nonpolys))


def msa_class_fit_chain(chain: Protein, prot_class: ProteinClass, custom_seq=None):
    prot_class_msa = msa_dict(
        phylogenetic_correction_taxid=chain.src_organism_ids[0],
        include_only_classes=[prot_class],
    )[prot_class]
    extended_msa = msaclass_extend_temp(
        prot_class,
        prot_class_msa,
        chain.entity_poly_seq_one_letter_code_can if custom_seq == None else custom_seq,
        chain.auth_asym_id,
        chain.parent_rcsb_id,
    )
    for seq in extended_msa:
        if prot_class in seq.getLabel():
            return str(seq)
    raise AssertionError(
        "Could not find sequence in for {} in aligned class. Likely label error {}".format(
            prot_class
        )
    )


def exit_port_posn(rcsb_id: str) -> list[float]:
    def __ul23():
        aln_conserved = [340, 351]
        return aln_conserved
        rcsb_id = "5aka"
        prot_class = "uL23"

        RA = RibosomeAssets(rcsb_id)
        chain = RA.get_prot_by_nomclass(prot_class)

        chainseq_ = RA.biopython_chain_get_seq(
            RA.biopython_structure(), chain.auth_asym_id, "protein"
        )
        print(chainseq_)

        prot_class_msa = msa_dict(
            phylogenetic_correction_taxid=chain.src_organism_ids[0],
            include_only_classes=[prot_class],
        )[prot_class]
        extended_msa = msaclass_extend_temp(
            prot_class,
            prot_class_msa,
            chainseq_,
            chain.auth_asym_id,
            chain.parent_rcsb_id,
        )

        for seq in extended_msa:
            if prot_class in seq.getLabel():
                print(seq.getLabel())
                subseq_ids = [74, 79]
                fwd_resids = [
                    util__forwards_match(str(seq), aligned_id)
                    for aligned_id in subseq_ids
                ]
                print(fwd_resids)
                print(SeqMatch.hl_ixs(str(seq), fwd_resids))
            else:
                print(seq.getLabel())
                # print(str(seq))
                print(SeqMatch.hl_ixs(str(seq), aln_conserved, color=92))
        return aln_conserved

    def __ul24():
        aln_conserved = [123, 124, 125]
        return aln_conserved
        rcsb_id = "5AKA"
        prot_class = "uL24"

        RA = RibosomeAssets(rcsb_id)
        chain = RA.get_prot_by_nomclass(prot_class)

        chainseq_ = RA.biopython_chain_get_seq(
            RA.biopython_structure(), chain.auth_asym_id, "protein"
        )
        print(chainseq_)

        prot_class_msa = msa_dict(
            phylogenetic_correction_taxid=chain.src_organism_ids[0],
            include_only_classes=[prot_class],
        )[prot_class]
        extended_msa = msaclass_extend_temp(
            prot_class,
            prot_class_msa,
            chainseq_,
            chain.auth_asym_id,
            chain.parent_rcsb_id,
        )

        for seq in extended_msa:
            if prot_class in seq.getLabel():
                print(seq.getLabel())
                subseq_ids = [46, 47, 48]
                fwd_resids = [
                    util__forwards_match(str(seq), aligned_id)
                    for aligned_id in subseq_ids
                ]
                print(SeqMatch.hl_ixs(str(seq), fwd_resids))
            else:
                print(seq.getLabel())
                print(SeqMatch.hl_ixs(str(seq), aln_conserved, color=92))

        return aln_conserved

    ra = RibosomeAssets(rcsb_id)
    bp_struct = ra.biopython_structure()

    ul23 = ra.get_prot_by_nomclass("uL23")
    ul24 = ra.get_prot_by_nomclass("uL24")

    if ul23 == None or ul24 == None:
        raise LookupError("Could not find uL23 or uL24 in {}".format(rcsb_id))

    print(ul23.auth_asym_id)
    print(ul24.auth_asym_id)

    bp_ul23_seq = ra.biopython_chain_get_seq(bp_struct, ul23.auth_asym_id, "protein")
    bp_ul24_seq = ra.biopython_chain_get_seq(bp_struct, ul24.auth_asym_id, "protein")

    bp_ul23 = ra.biopython_get_chain(ul23.auth_asym_id)
    bp_ul24 = ra.biopython_get_chain(ul24.auth_asym_id)

    aligned_ul23 = msa_class_fit_chain(ul23, "uL23", custom_seq=bp_ul23_seq)
    aligned_ul24 = msa_class_fit_chain(ul24, "uL24", custom_seq=bp_ul24_seq)

    backwards_mapped_ul24 = [
        util__backwards_match(aligned_ul24, residue) for residue in __ul24()
    ]
    backwards_mapped_ul23 = [
        util__backwards_match(aligned_ul23, residue) for residue in __ul23()
    ]

    residues_ul23 = [bp_ul24[i] for i in backwards_mapped_ul24]
    residues_ul24 = [bp_ul23[i] for i in backwards_mapped_ul23]

    ul23_M = np.average([*map(lambda res: res.center_of_mass(), residues_ul23)], axis=0)
    ul24_M = np.average([*map(lambda res: res.center_of_mass(), residues_ul24)], axis=0)

    return np.average([ul23_M, ul24_M], axis=0).tolist()


def make_cylinder(p1: list[float], p2: list[float], R: float):
    height = math.sqrt(
        (p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2 + (p2[2] - p1[2]) ** 2
    )
    center = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)

    return {"center": center, "height": height, "radius": R}


def pt_is_inside_cylinder(cylinder, point):
    distance_xy = math.sqrt(
        (point[0] - cylinder["center"][0]) ** 2
        + (point[1] - cylinder["center"][1]) ** 2
    )

    if distance_xy <= cylinder["radius"]:
        distance_top = abs(point[2] - (cylinder["center"][2] + cylinder["height"] / 2))
        distance_bottom = abs(
            point[2] - (cylinder["center"][2] - cylinder["height"] / 2)
        )

        if (
            distance_top <= cylinder["height"] / 2
            and distance_bottom <= cylinder["height"] / 2
        ):
            return True

    return False


def _x_euk_structs():
    EUK_STRUCTS = []
    with open("eukarya_03_07_2023.txt", "r") as data_file:
        for line in data_file:
            structs = line.split(",")
            EUK_STRUCTS = [*EUK_STRUCTS, *structs]
    return EUK_STRUCTS


def _x_obtain_euk_structs():
    EUK_STRUCTS = []
    with open("eukarya_2023.txt", "r") as data_file:
        for line in data_file:
            structs = line.split()
            EUK_STRUCTS = [*EUK_STRUCTS, *structs]

    RCSB_ID = "4UG0"

    al = Assetlist(
        profile                 = True,
        structure               = True,
        structure_modified      = False,
        chains_and_modified_cif = False,
        factors_and_ligands     = False,
        png_thumbnail           = False,
    )

    obtain_assets_threadpool(EUK_STRUCTS, al)


if __name__ == "__main__":
    import numpy as np

    EUK        = _x_euk_structs()
    PTC_COORDS = {}
    for RCSB_ID in EUK :
        try:
            print("---------+--------------")
            print("Processing {}".format(RCSB_ID))
            ress, auth_asym_id = ptc_resdiues_get(RCSB_ID)
            midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)

            residue_labels = [(res.get_resname(), res.id[1]) for res in ress]
            print(residue_labels)

            writeout = {
                "site_9_residues": [
                    (res.get_resname(), res.get_segid()) for res in ress
                ],
                "LSU_rRNA_auth_asym_id": auth_asym_id,
                "midpoint_coordinates": midpoint_coords,
            }

            PTC_COORDS = {**PTC_COORDS, RCSB_ID: writeout}

        except Exception as e:
            print(e)

    # with open("ptc_coords_eukarya.json", 'w') as outfile:
    #     json.dump(PTC_COORDS, outfile, indent=4)

    # asyncio.run( obtain_assets(RCSB_ID, al))

    # for rcsb_id in EUKARYOTIC:
    # - identify the 25/28/35S
    # - fuzzyfind the conserved residues
    # - triangulate from the tips of the [conserved res.] comb

# Edge cases:
# 8i9x is a preribosome with too many gaps to identify ptc comb from doris-site-9 (should try with more)