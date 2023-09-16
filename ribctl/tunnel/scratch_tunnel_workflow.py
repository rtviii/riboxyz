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
