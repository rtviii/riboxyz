import math
import os
from fuzzysearch import find_near_matches
import numpy as np
from ribctl.lib.schema.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
    ResidueSummary,
)
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from functools import reduce
from Bio.PDB.Structure import Structure
from ribctl.lib.schema.types_ribosome import RNA
from fuzzysearch import find_near_matches


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


def ptc_fuzzyfind_subseq_in_chain( biopython_struct, auth_asym_id: str, assembly_id: int = 0 )  \
    -> tuple[list[Residue], list[Residue], list[Residue], str]:

    chain3d       : Chain         = biopython_struct.child_dict[assembly_id].child_dict[auth_asym_id]
    ress          : list[Residue] = chain3d.child_list
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

    match9 = pick_match(
        find_near_matches(DORIS_ET_AL["SITE_9"], raw_seq, max_l_dist=1), len(raw_seq)
    )
    match8 = pick_match(
        find_near_matches(DORIS_ET_AL["SITE_8"], raw_seq, max_l_dist=1), len(raw_seq)
    )
    match6 = pick_match(
        find_near_matches(DORIS_ET_AL["SITE_6"], raw_seq, max_l_dist=1), len(raw_seq)
    )

    PTC_residues_9 = (
        [ress_sanitized[i] for i in list(range(match9.start, match9.end))]
        if match9
        else []
    )
    PTC_residues_8 = (
        [ress_sanitized[i] for i in list(range(match8.start, match8.end))]
        if match8
        else []
    )
    PTC_residues_6 = (
        [ress_sanitized[i] for i in list(range(match6.start, match6.end))]
        if match6
        else []
    )

    return PTC_residues_6, PTC_residues_8, PTC_residues_9, auth_asym_id


def ptc_resdiues_get(
    biopython_structure: Structure, rnas: list[RNA], assembly_id: int = 0
) -> tuple[list[Residue], str]:
    """
    Given a bacterial(?) structure (assertion un_enforced):
     - obtain the LSU rRNA chain
     - sanitze the sequence (canonicalize non-canonical residues)
     - locate the conserved subsequence (due to Doris et al. 2015)
     Returns the residue list for the ptc and the `auth_asym_id` of the rRNA chain
    """
    struct_profile = biopython_structure
    # TODO: This whole mess has to go. The LSU_rRNA will be much easier to identify once we have the rRNA HMHs in place.
    # For now just know that this is a hack.

    matches = {}
    for p_rna in rnas:
        m6, m8, m9, auth_asym_id = ptc_fuzzyfind_subseq_in_chain( struct_profile, p_rna.auth_asym_id )
        matches = {**matches, auth_asym_id: [m6, m8, m9]}

    try:
        auth_asym_id, rRNA_fragment_matches = list( filter(lambda match_kv: len(match_kv[1][2]) > 0, list(matches.items())) )[0]
    except Exception as e:
        raise Exception("Error:Could not identify PTC residues in {} : {} ".format(biopython_structure.id, e)   )

    chain3d: Chain = struct_profile.child_dict[assembly_id].child_dict[auth_asym_id]
    ress: list[Residue] = chain3d.child_list

    ress_sanitized: list[Residue] = [
        *filter(lambda r: r.get_resname() in ["A", "C", "G", "U", "PSU"], ress)
    ]

    for _r in ress_sanitized:
        if _r.get_resname() == "PSU":
            _r.resname = "U"

    raw_seq      = reduce(lambda x, y: x + y.resname, ress_sanitized, "")

    matches      = find_near_matches(DORIS_ET_AL["SITE_9"], raw_seq, max_l_dist=1)
    m0           = pick_match(matches, len(raw_seq))
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
