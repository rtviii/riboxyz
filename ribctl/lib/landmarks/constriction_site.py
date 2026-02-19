# ribctl/lib/landmarks/constriction_site.py

from __future__ import annotations

from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Residue import Residue
import numpy as np
from loguru import logger

from ribctl.lib.exceptions import SkipAsset
from ribctl.lib.types.polymer.base import CytosolicProteinClass, MitochondrialProteinClass
from ribctl.lib.utils import find_closest_pair_two_sets, midpoint
from ribctl.ribosome_ops import RibosomeOps


def _best_matching_protein(profile, accepted_classes: set) -> object | None:
    """
    Return the best matching protein Polymer from profile.proteins whose nomenclature
    contains any class in accepted_classes. Prefer longest sequence as a tiebreaker.
    """
    proteins = getattr(profile, "proteins", None) or []
    candidates = []
    for p in proteins:
        try:
            noms = set(p.nomenclature or [])
        except Exception:
            continue
        if any(n in accepted_classes for n in noms):
            candidates.append(p)

    if not candidates:
        return None

    # Prefer longest sequence; fallback to first if any missing lengths
    candidates.sort(key=lambda x: getattr(x, "entity_poly_seq_length", 0), reverse=True)
    return candidates[0]


def _get_chain_any_model(structure: Structure, auth_asym_id: str) -> Chain:
    """
    Try to retrieve chain auth_asym_id from any model in the Bio.PDB structure.
    """
    for model in structure:
        if isinstance(model, Model) and auth_asym_id in model.child_dict:
            return model[auth_asym_id]
    raise SkipAsset(f"Chain {auth_asym_id} not found in any model of mmCIF structure")


def _safe_center_of_mass(res: Residue) -> np.ndarray | None:
    """
    Compute residue center of mass; return None if residue has no atoms or computation fails.
    """
    try:
        if len(res.child_list) == 0:
            return None
        return res.center_of_mass()
    except Exception:
        return None


def get_constriction(rcsb_id: str) -> np.ndarray:
    ro = RibosomeOps(rcsb_id)
    profile = ro.profile

    # Determine which nomenclature keys to seek
    if profile.mitochondrial:
        ul4_keys = {MitochondrialProteinClass.uL4m}
        ul22_keys = {MitochondrialProteinClass.uL22m}
    else:
        # In eukaryotes, "uL22" might appear as "eL22" in some annotations.
        ul4_keys = {CytosolicProteinClass.uL4}
        ul22_keys = {CytosolicProteinClass.uL22, CytosolicProteinClass.eL22}

    uL4 = _best_matching_protein(profile, ul4_keys)
    uL22 = _best_matching_protein(profile, ul22_keys)

    if uL4 is None or uL22 is None:
        raise SkipAsset(f"Could not find uL4/uL22 (or eL22) in {rcsb_id} — likely SSU-only or incomplete LSU")

    structure = ro.assets.biopython_structure()

    # Resolve chains robustly across models
    try:
        uL4_c: Chain = _get_chain_any_model(structure, uL4.auth_asym_id)
        uL22_c: Chain = _get_chain_any_model(structure, uL22.auth_asym_id)
    except KeyError as e:
        raise SkipAsset(f"Chain not found while locating constriction site in {rcsb_id}: {e}")

    uL4_coords = []
    for r in uL4_c.child_list:
        com = _safe_center_of_mass(r)
        if com is not None:
            uL4_coords.append(com)

    uL22_coords = []
    for r in uL22_c.child_list:
        com = _safe_center_of_mass(r)
        if com is not None:
            uL22_coords.append(com)

    if len(uL4_coords) == 0 or len(uL22_coords) == 0:
        raise SkipAsset(f"Insufficient coordinates for uL4/uL22 in {rcsb_id} (empty residue COM list)")

    try:
        p1, p2 = find_closest_pair_two_sets(uL4_coords, uL22_coords)
    except Exception as e:
        raise SkipAsset(f"Failed closest-pair search for constriction site in {rcsb_id}: {e}")

    return midpoint(p1, p2)

