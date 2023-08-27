import json
import math
from ribctl.lib.msalib import AMINO_ACIDS_3_TO_1_CODE, msa_dict, msaclass_extend_temp, util__backwards_match, util__forwards_match
from ribctl import RIBETL_DATA
from ribctl.lib.types.types_ribosome import PolymericFactor, Protein, ProteinClass
from ribctl.lib.types.types_binding_site import AMINO_ACIDS, NUCLEOTIDES, ResidueSummary
from ribctl.ribosome_assets import RibosomeAssets
import loguru
from Bio.PDB.Atom import Atom
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from functools import reduce
import numpy as np
import inspect
import os
from colorama import Fore, init

from api.ribctl.lib.mod_transpose_bsites import SeqMatch
init(autoreset=True)

BACTERIAL = ['4V4I', '4WSD', '4V87', '4WQ1', '4V5D', '1NWX', '4V6Z', '6QNR', '5T7V', '4V9J', '5GAE', '4V4Q', '7UVZ', '7ZQ6', '7QV3', 
              '5LZF', '6GXO', '5MYJ', '4V55', '7OF4', '6HMA', '7UW1', '7RQD', '7K53', '6GSL', '4W2G', '6UO1', '7MSM', '4WRO', '4LEL',
              '5L3P', '4V48', '7N1P', '4V6N', '4WQU', '4WFN', '7K51', '5O61', '8B7Y', '7P3K', '6ORL', '4V6Y', '6XDQ', '6Q9A', '5EL4',
              '6S0K', '6WDD', '5U4I', '7JSW', '3JCN', '4V4Z', '4WU1', '7RQA', '6SPG', '1VVJ', '7MSZ', '4V6Q', '4V9Q', '4V72', '4V75',
              '6GSK', '5KCS', '7UPH', '6V3D', '4V8C', '5WE6', '4V4Y', '6SPD', '7MT2', '6WNT', '5HCP', '4WQY', '4V5H', '1ML5', '5EL6',
              '4WSM', '6OF6', '7PJX', '6XQD', '5DM6', '5VP2', '8A57', '4V4N', '4V6F', '4V7M', '3JAH', '8B0X', '5KPS', '4IOC', '7JT3',
              '8AYE', '4V97', '3J7Z', '6WD3', '5LZB', '4V9P', '8EKB', '6QNQ', '6ENJ', '5ADY', '4V5Q', '5J88', '4V42', '1NKW', '4LT8',
              '5ND9', '5LZW', '5DM7', '4W29', '3J5L', '5DFE', '4V6O', '4V5S', '5M1J', '4V51', '4L71', '1VY5', '7K00', '3J9Y', '5LZD',
              '7ZQ5', '7UVY', '6OM6', '7A5G', '5A9Z', '7OJ0', '8A63', '6WDK', '6GXN', '8EIU', '7SS9', '7N2V', '6BU8', '7JT2', '5CZP',
              '4XEJ', '6FKR', '4V95', '6BOH', '6VYR', '6ND6', '3JBU', '3J9Z', '4V4A', '5ZET', '4U1U', '6WD7', '6NDK', '6GXM', '4V4H',
              '7ZOD', '7UVX', '5AFI', '4V5F', '6Y69', '4P6F', '6V3A', '7JSS', '4YZV', '7OII', '7NHL', '6QDW', '6HA8', '4V8T', '5WFK',
              '6V3B', '4V5Y', '3JCJ', '4V8U', '5IBB', '6OFX', '4IO9', '5LZZ', '7PAU', '7P6Z', '7MSC', '5AKA', '4V7L', '4V5E', '4WQF',
              '7N2C', '5HL7', '4V5O', '4W2I', '7P7U', '6CFK', '4V5R', '7YLA', '7RQC', '4Y4O', '3BBX', '4V63', '6N1D', '6HTQ', '4V9O',
              '3DLL', '5UQ7', '5H5U', '4WZO', '4V54', '6O8Y', '6BUW', '6NSH', '4V8O', '6BY1', '6BZ6', '6VZ3', '7M4W', '7ASM', '5UYL',
              '5WFS', '1NWY', '6WDA', '4V68', '5LZC', '6WD5', '6B4V', '4WZD', '6X7F', '7AZS', '4V7A', '4V7S', '6S0Z', '4UY8', '3JAI',
              '7A5F', '4V6T', '4V8J', '7MSH', '4V6D', '5IMQ', '5GAF', '5ZEP', '4V8A', '7KGB', '4YPB', '5VYC', '6TMF', '1YIT', '7P7Q',
              '7JIL', '7N2U', '5LI0', '4V83', '4V4X', '4TUE', '4WT1', '4V7W', '5LZX', '6WDM', '7SSN', '6WDL', '7ST6', '5GAG', '4V8B',
              '6YS3', '7MD7', '2ZJR', '7UG7', '4LSK', '6VU3', '4V65', '4V70', '7JQC', '6WOO', '6O97', '4WCE', '4V7J', '4V8H', '6OPE',
              '6WDC', '6HRM', '7NHN', '6CFJ', '6OGF', '6SPB', '4U26', '7RQB', '7U2J', '5MDV', '4V6G', '6DZI', '4V7I', '2RDO', '7K55',
              '4V5A', '7Q4K', '4V7X', '6BZ7', '7UNW', '5UYN', '7SSO', '7BV8', '8G61', '5U9F', '4V9I', '5WIS', '4WRA', '6GBZ', '4TUB',
              '6WD4', '5MDZ', '6ORE', '6OSK', '4V9N', '4V4G', '7NHK', '4V90', '5LZT', '4LFZ', '4V6L', '4U1V', '5WF0', '5GAD', '6WNW',
              '4V89', '7OIG', '5EL7', '4TUD', '6U48', '4V9B', '6VYQ', '7QGH', '6OT3', '6OTR', '4W2F', '7AQD', '8CVJ', '6VYW', '4CSU',
              '5J30', '4V7V', '6WD0', '4Y4P', '6Z6K', '5EL5', '5U9G', '4U25', '5KCR', '6VYZ', '7LH5', '6DZP', '6OSQ', '7M5D', '5KPX',
              '7PJY', '7RYG', '6S13', '1VY6', '5AA0', '6X9Q', '7NWT', '4V6E', '6VZ5', '1W2B', '6WDI', '7PJZ', '7OT5', '6OGI', '8CVL',
              '5J3C', '7M4V', '7RYF', '7OOD', '6I0Y', '7JQB', '5ZLU', '5MDY', '7SFR', '5XYM', '7OIF', '4V84', '6VWL', '6VZ2', '7N30',
              '6WD9', '6SZS', '4V6V', '5NDK', '8CVK', '7BL4', '4V57', '5APO', '7PJW', '4WOI', '5J8B', '4L47', '6ENF', '4V8F', '6O90',
              '6XHX', '7P7R', '7O5B', '3CF5', '7S1J', '7ZTA', '6XZA', '6WDE', '6OST', '6O9J', '6WD8', '4U27', '4V8E', '6VYU', '6WD1',
              '7JSZ', '5KPV', '4TUC', '6OSI', '5J4D', '7QGN', '5E7K', '7RQ8', '6OUO', '2ZJP', '4V79', '6TC3', '4V5K', '4V66', '4WFB',
              '4V4T', '3PIO', '5IMR', '5V7Q', '4V7K', '4V8D', '8G6Y', '6YSU', '7SA4', '5DOY', '4V7Y', '4U24', '4V5J', '5MGP', '4WWW',
              '6CFL', '4V78', '8FOM', '7PJT', '4V6A', '6GZQ', '5E81', '6BOK', '7ASO', '4V5P', '5UYM', '4V77', '4U67', '7D6Z', '4V71',
              '5JVG', '5JVH', '4V64', '7SSL', '6WDH', '4WT8', '4V6S', '4V8X', '4V52', '7U2H', '6VZJ', '7QG8', '5LZY', '7P7S', '7M4Y',
              '6XDR', '7LV0', '7BL5', '6ENU', '6XHV', '4V8Q', '7NSO', '4Z3S', '5APN', '6XQE', '1XBP', '5J4B', '7S1G', '5UQ8', '6VYS',
              '5IB8', '4V4V', '6XZ7', '6N9E', '6TBV', '7AZO', '6XZB', '6O9K', '5JTE', '7AQC', '4V6K', '4V9S', '6S12', '5NCO', '6GC0',
              '1YJW', '4V7T', '7PJU', '6WRU', '7UVW', '6ORD', '7OPE', '4V5B', '4V7Z', '6W6P', '7S1K', '6OG7', '6YSS', '6O8Z', '5GAH',
              '4V74', '4WF1', '4V49', '6O3M', '5UYP', '4P70', '5W4K', '7QGU', '4V4W', '6H58', '6I7V', '4V5C', '1SM1', '6VYT', '4W4G',
              '4WR6', '4U20', '7B5K', '6CAE', '6GWT', '7PJS', '5VPO', '6YSI', '5UYK', '5NWY', '4V69', '7S1H', '5J4C', '4IOA', '7BL2',
              '6XHY', '7UNU', '3JCE', '5LZA', '7OIZ', '7NSQ', '3JAG', '7SSD', '4V85', '6XHW', '6VYY', '2J28', '6WNV', '7RYH', '6DNC',
              '4V8G', '6VWM', '6WDG', '6VYX', '4V8I', '7M4Z', '4V5G', '4BTS', '6C5L', '7TOS', '4WQR', '5IQR', '6X7K', '6GZZ', '6YEF',
              '7NWW', '5JC9', '6H4N', '7P7T', '7BL3', '3JCD', '6O8W', '5V8I', '5WE4', '6O8X', '7UNV', '5WDT', '4LNT', '4V47', '6WD2',
              '3J9W', '8A5I', '7K50', '4V6P', '6GXP', '4V9C', '7S1I', '4V4J', '4V9R', '4YBB', '6OJ2', '7F0D', '6NUO', '7JT1', '7NHM',
              '6WDF', '7ASP', '4V4P', '4V4R', '7ST2', '6S0X', '5IT8', '7S0S', '7OF2', '6NTA', '4V9H', '7U2I', '4V9D', '8G6X', '7UNR',
              '5NGM', '6VZ7', '4V5L', '6Q98', '3JA1', '8G6W', '4V7U', '6Q97', '8EKC', '4W2E', '6YSR', '4V67', '6CZR', '5O60', '4V5N',
              '5WIT', '5OT7', '5J7L', '7QH4', '4WFA', '7UVV', '5FDV', '7OTC', '3J8G', '5J8A', '4V6C', '5LZV', '6WDB', '4V9A', '5JU8',
              '6V39', '7QQ3', '7LVK', '7BL6', '4V4S', '6C4I', '4V7B', '6UCQ', '4V53', '5FDU', '4V9K', '7RQ9', '5VPP', '4ZSN', '7K52',
              '7NSP', '8FON', '5DOX', '5UYQ', '4V50', '6ND5', '6VWN', '3JBV', '6HA1', '6OGG', '7QV2', '8BUU', '4V9L', '7PAT', '5ZEB',
              '5NDJ', '5TCU', '7QGR', '6PJ6', '7MT3', '7RQE', '6GC8', '6Q95', '7N31', '6OXA', '4V9M', '6OF1', '5J91', '6G5I', '7D80',
              '2ZJQ', '4V5M', '4V6R', '7QV1', '6N9F', '4V73', '4WPO', '6WD6', '1VY4', '6YHS', '6WDJ', '5LZE', '6GZX', '5J5B', '7K54',
              '8G5Z', '6FXC', '6NWY', '6YST', '5V93', '4V56', '8C8X', '7MT7', '7PJV', '7ST7', '5KPW', '4W2H', '6SPF', '7A5J', '4V7P',
              '4V8N', '7M4X', '4V76', '5MDW', '4TUA', '3PIP', '6GSJ', '5ND8', '7SSW', '6BZ8', '1VY7']

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4574749/pdf/1719.pdf
DORIS_ET_AL = {
    "subseq_6": "AAGACCC",
    "subseq_8": "GGAUAAC",
    "subseq_9": "GAGCUGGGUUUA",
    'b': {
        "site_6": [2059, 2060, 2061, 2062, 2063, 2064, 2065],
        "site_8": [2446, 2447, 2448, 2449, 2450, 2451, 2452],
        "site_9": [2576, 2577, 2578, 2579, 2580, 2581, 2582, 2583, 2584, 2585, 2586, 2587],
    }
}

def pick_match(matches, rna_length: int):
    if len(matches) == 0:
        print("No matches found!!")
        exit(1)

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

def residue_is_canonical(res: Residue | ResidueSummary) -> bool:
    return res.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]

def ptc_resdiues_get(rcsb_id: str, assembly_id: int = 0) -> tuple[list[Residue], str]:
    """Given a bacterial structure (assertion unenforced):
     - obtain the LSU rRNA chain
     - sanitze the sequence (canonicalize non-canonical residues)
     - locate the conserved subsequence (due to Doris et al.)
     - return the residue list for the ptc
     """

    R = RibosomeAssets(rcsb_id)
    struct_profile = R.biopython_structure()
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

    from fuzzysearch import find_near_matches

    matches = find_near_matches(DORIS_ET_AL["subseq_9"], raw_seq, max_l_dist=1)
    m0 = pick_match(matches, len(raw_seq))
    PTC_residues = [ress_sanitized[i] for i in list(range(m0.start, m0.end))]

    return PTC_residues, auth_asym_id

def ptc_residues_to_atom_coordinates(
    reslist: list[Residue],
    auth_asym_id: str
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
                map(lambda x: float(x), list(atom_coords)))

    ptc_coordinates = {auth_asym_id: ptc_coordinates}

    return ptc_coordinates

def ptc_residues_calculate_midpoint(reslist: list[Residue], auth_asym_id: str) -> list[float]:
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
    subseq_length   = len(subseq_residues)

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

def tunnel_obstructions(rcsb_id: str, ptc_midpoint: tuple[float, float, float], radius: int = 30) -> tuple[list[PolymericFactor], list[ResidueSummary]]:
    R = RibosomeAssets(rcsb_id)
    structure, profile = R.get_struct_and_profile()

    neigbor_search = NeighborSearch(list(structure.get_atoms()))
    nbr_residues = []

    nbr_residues.extend(neigbor_search.search(ptc_midpoint, radius, level='R'))
    nbr_residues = list(
        set([* map(ResidueSummary.from_biopython_residue, nbr_residues)]))

    foreign_nonpolys = []
    foregin_polymers = []

    for N in nbr_residues:

        if not residue_is_canonical(N):
            foreign_nonpolys.append(N)
        parent_chain, chain_type = R.get_chain_by_auth_asym_id(
            N.parent_auth_asym_id)
        if (parent_chain, chain_type) != (None, None):
            if chain_type == "PolymericFactor":
                foregin_polymers.append(parent_chain)
        else:
            raise LookupError(
                "Parent chain must exist in structure. Something went wrong (with the data, probably)")

    print("Inspected midpoint ", ptc_midpoint, " with radius", radius)
    return list(set(foregin_polymers)), list(set(foreign_nonpolys))

def msa_class_fit_chain(chain: Protein, prot_class: ProteinClass, custom_seq=None):
    prot_class_msa = msa_dict(phylogenetic_correction_taxid=chain.src_organism_ids[0], include_only_classes=[prot_class])[prot_class]
    extended_msa   = msaclass_extend_temp(prot_class, prot_class_msa, chain.entity_poly_seq_one_letter_code_can if custom_seq == None else custom_seq, chain.auth_asym_id, chain.parent_rcsb_id)
    for seq in extended_msa:
        if prot_class in seq.getLabel():
            return str(seq)
    raise AssertionError(
        "Could not find sequence in for {} in aligned class. Likely label error {}".format(prot_class))

def exit_port_posn(rcsb_id:str)->list[float]:

    def __ul23():
        aln_conserved = [340, 351]
        return aln_conserved
        rcsb_id    = '5aka'
        prot_class = 'uL23'

        RA = RibosomeAssets(rcsb_id)
        chain = RA.get_prot_by_nomclass(prot_class)

        chainseq_ = RA.biopython_chain_get_seq(
            RA.biopython_structure(), chain.auth_asym_id, 'protein')
        print(chainseq_)

        prot_class_msa = msa_dict(
            phylogenetic_correction_taxid=chain.src_organism_ids[0], include_only_classes=[prot_class])[prot_class]
        extended_msa = msaclass_extend_temp(
            prot_class, prot_class_msa, chainseq_, chain.auth_asym_id, chain.parent_rcsb_id)

        for seq in extended_msa:
            if prot_class in seq.getLabel():
                print(seq.getLabel())
                subseq_ids = [74, 79]
                fwd_resids = [util__forwards_match(
                    str(seq), aligned_id) for aligned_id in subseq_ids]
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
        rcsb_id = '5AKA'
        prot_class = 'uL24'

        RA = RibosomeAssets(rcsb_id)
        chain = RA.get_prot_by_nomclass(prot_class)

        chainseq_ = RA.biopython_chain_get_seq(
            RA.biopython_structure(), chain.auth_asym_id, 'protein')
        print(chainseq_)

        prot_class_msa = msa_dict(
            phylogenetic_correction_taxid=chain.src_organism_ids[0], include_only_classes=[prot_class])[prot_class]
        extended_msa = msaclass_extend_temp(
            prot_class, prot_class_msa, chainseq_, chain.auth_asym_id, chain.parent_rcsb_id)

        for seq in extended_msa:
            if prot_class in seq.getLabel():
                print(seq.getLabel())
                subseq_ids = [46, 47, 48]
                fwd_resids = [util__forwards_match(
                    str(seq), aligned_id) for aligned_id in subseq_ids]
                print(SeqMatch.hl_ixs(str(seq), fwd_resids))
            else:
                print(seq.getLabel())
                print(SeqMatch.hl_ixs(str(seq), aln_conserved, color=92))

        return aln_conserved
    ra        = RibosomeAssets(rcsb_id)
    bp_struct = ra.biopython_structure()

    ul23      = ra.get_prot_by_nomclass('uL23')
    ul24      = ra.get_prot_by_nomclass('uL24')

    if ul23 == None or ul24 == None:
        raise LookupError("Could not find uL23 or uL24 in {}".format(rcsb_id))

    print(ul23.auth_asym_id)
    print(ul24.auth_asym_id)

    bp_ul23_seq = ra.biopython_chain_get_seq(bp_struct, ul23.auth_asym_id, 'protein')
    bp_ul24_seq = ra.biopython_chain_get_seq(bp_struct, ul24.auth_asym_id, 'protein')
    
    bp_ul23 = ra.biopython_get_chain(ul23.auth_asym_id)
    bp_ul24 = ra.biopython_get_chain(ul24.auth_asym_id)

    aligned_ul23 = msa_class_fit_chain(ul23, 'uL23', custom_seq=bp_ul23_seq)
    aligned_ul24 = msa_class_fit_chain(ul24, 'uL24', custom_seq=bp_ul24_seq)

    backwards_mapped_ul24 = [util__backwards_match(aligned_ul24, residue) for residue in __ul24()]
    backwards_mapped_ul23 = [util__backwards_match(aligned_ul23, residue) for residue in __ul23()]

    residues_ul23 = [bp_ul24[i] for i in backwards_mapped_ul24]
    residues_ul24 = [bp_ul23[i] for i in backwards_mapped_ul23]

    ul23_M = np.average([*map(lambda res:res.center_of_mass(), residues_ul23)], axis=0)
    ul24_M = np.average([*map(lambda res:res.center_of_mass(), residues_ul24)], axis=0)

    return np.average([ul23_M, ul24_M],axis=0).tolist()

def make_cylinder(p1:list[float], p2:list[float], R:float):
    height = math.sqrt((p2[0] - p1[0])**2 + (p2[1] - p1[1])**2 + (p2[2] - p1[2])**2)
    center = ((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2)
    
    return {
        "center": center,
        "height": height,
        "radius": R
    }

def pt_is_inside_cylinder(cylinder, point):

    distance_xy = math.sqrt((point[0] - cylinder["center"][0])**2 + (point[1] - cylinder["center"][1])**2)

    if distance_xy <= cylinder["radius"]:

        distance_top    = abs(point[2] - (cylinder["center"][2] + cylinder["height"]/2))
        distance_bottom = abs(point[2] - (cylinder["center"][2] - cylinder["height"]/2))

        if distance_top <= cylinder["height"]/2 and distance_bottom <= cylinder["height"]/2:
            return True

    return False

if __name__ == "__main__":

    import numpy as np
    import matplotlib.pyplot as plt


    # rcsb_id = "4V4I"
    # struct,profile = RibosomeAssets(rcsb_id).get_struct_and_profile()
    # atoms  = struct.get_atoms()
    # cylinder = make_cylinder([-170.46199798583984, 20.90049934387207, -82.8115005493164], [-192.08737182617188, 70.13729095458984, -36.43252182006836], 7)
    # # tunnel_atoms = [*filter(lambda atom: point_inside_cylinder(cylinder, atom.get_coord()), atoms)]
    # for atom in atoms:
    #     print(atom.get_coord())
    #     print(point_inside_cylinder(cylinder, list(atom.get_coord())))
    # print(len(atoms))
    # print(len(tunnel_atoms))

    # for rcsb_id in BACTERIAL:
    #     try:
    #         PTC_COMB_PATH  = os.path.join("/home/rxz/dev/docker_ribxz/api/ribctl/assets/landmarks", "{}.json".format(rcsb_id))
    #         if os.path.exists(PTC_COMB_PATH):
    #             print("File exists, skipping {}".format(PTC_COMB_PATH))
    #             continue
    #         posn_exit_port = exit_port_posn(rcsb_id)
    #         posn_ptc       = ptc_residues_calculate_midpoint(*ptc_resdiues_get(rcsb_id, 0))

    #         landmarks_dict = {
    #             "ptc"      : posn_ptc,
    #             "exit_port": posn_exit_port
    #         }


    #         with open(PTC_COMB_PATH, 'w') as f:
    #             json.dump(landmarks_dict, f)
    #             print("Saved {}".format(PTC_COMB_PATH))


    #     except Exception as e:
    #         print("Failed to process {}".format(rcsb_id))
    #         print(e)
    # ?---------------------------------------------------?#

    # msaclass_extend_temp()

    # for RCSB_ID in BACTERIAL:
    #     print("Processing {}".format(RCSB_ID))

    # try:
    #     PTC_COMB_PATH = os.path.join("/home/rxz/dev/docker_ribxz/api/ribctl/assets/ptc_comb_midpoint", "{}.json".format(RCSB_ID))
    #     residues,auth_asym_id = ptc_resdiues_get(RCSB_ID, 0)
    #     ptc_coordinates = ptc_residues_to_atom_coordinates(residues, auth_asym_id)
    #     print("Saving {}".format(PTC_COMB_PATH))

    # except Exception as e:
    #     loguru.logger.error("Failed to extract PTC {}.\n{}".format(RCSB_ID, e))
    #     ...

    # print("Writing to ", PTC_RESIDUES_PATH)

    # ptcres, auth_asym_id = ptc_resdiues_get(RCSB_ID, 0)
    # midpoint = ptc_midpoint(ptcres, auth_asym_id)
    # poly, nonpoly = tunnel_obstructions(RCSB_ID, midpoint)

    # print("-----------POLYMERICS-----------")
    # for pl in poly:
    #     pprint(pl)

    # print("-----------NONPOLYMERICS-----------")
    # for npl in nonpoly:
    #     pprint(npl)
