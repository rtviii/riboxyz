# ribctl/lib/npet2/core/structure_selection.py
"""
Chain selection policies for tunnel wall definition.

Works with npet2-internal RibosomeProfile types exclusively.
"""

from __future__ import annotations

from typing import Iterable, Set

from .ribosome_types import RibosomeProfile, PolymerEntry
from .polymer_enum import tRNA as tRNAClass

_TRNA_HINTS = ("trna", "transfer rna")


def looks_like_trna(poly: PolymerEntry) -> bool:
    """Heuristic tRNA detection from nomenclature and description."""
    for nom in poly.nomenclature:
        if isinstance(nom, tRNAClass):
            return True
        if "trna" in str(nom).lower():
            return True

    desc = (poly.rcsb_pdbx_description or "").lower()
    if any(h in desc for h in _TRNA_HINTS):
        return True

    return False


def ribosome_wall_auth_asym_ids(
    profile: RibosomeProfile,
    *,
    exclude_trna: bool = True,
    extra_exclude: Iterable[str] = (),
) -> Set[str]:
    """
    Tunnel wall = ribosomal proteins + rRNAs.

    This automatically excludes waters/ions/ligands/nonpolymers (they aren't
    in proteins or rnas). Modified residues within ribosomal polymers are
    included because they are covalently part of the wall.
    """
    wall = {p.auth_asym_id for p in profile.proteins}
    wall |= {r.auth_asym_id for r in profile.rnas}

    if exclude_trna:
        all_polys = profile.all_polymers()
        trna_ids = {p.auth_asym_id for p in all_polys if looks_like_trna(p)}
        wall -= trna_ids

    wall -= set(extra_exclude)
    return wall


def intersect_with_first_assembly(
    profile: RibosomeProfile, chain_ids: Set[str]
) -> Set[str]:
    """Only keep chains present in the first assembly (if assembly_map available)."""
    try:
        asm_ids = set(profile.first_assembly_auth_asym_ids())
        return chain_ids & asm_ids
    except (ValueError, IndexError):
        return chain_ids


def tunnel_debris_chains(rcsb_id: str, profile: RibosomeProfile) -> list[str]:
    """
    Hardcoded per-structure chain exclusions for known tunnel debris.

    Also handles mitochondrial mL45 (sits inside the exit tunnel).
    """
    from .polymer_enum import MitochondrialProteinClass

    DEBRIS_MAP = {
        "3J7Z": ["a", "7"],
        "5GAK": ["z"],
        "5NWY": ["s"],
        "7A5G": ["Y2"],
        "9F1D": ["BK"],
    }
    skip = list(DEBRIS_MAP.get(rcsb_id.upper(), []))

    if profile.mitochondrial:
        for poly in profile.all_polymers():
            if MitochondrialProteinClass.mL45 in poly.nomenclature:
                skip.append(poly.auth_asym_id)
                break

    return skip


def atom_inclusion_policy(
    profile: RibosomeProfile,
    config,
    rcsb_id: str,
) -> dict:
    """
    Central policy for which atoms go into occupancy calculations.

    Returns:
        wall_chain_ids: set of auth_asym_ids for occupancy
        excluded_chain_ids: set of auth_asym_ids removed
        reasons: dict mapping excluded chain_id -> reason string
    """
    debris = tunnel_debris_chains(rcsb_id, profile)
    manual_exclude = list(getattr(config, "occupancy_exclude_auth_asym_ids", ()))
    exclude_trna = bool(getattr(config, "occupancy_exclude_trna", True))

    all_exclude = list(dict.fromkeys(debris + manual_exclude))

    wall = ribosome_wall_auth_asym_ids(
        profile, exclude_trna=exclude_trna, extra_exclude=all_exclude
    )
    wall = intersect_with_first_assembly(profile, wall)

    reasons = {}
    for c in debris:
        reasons[c] = "tunnel_debris (hardcoded)"
    for c in manual_exclude:
        if c not in reasons:
            reasons[c] = "config exclude"
    if exclude_trna:
        for p in profile.all_polymers():
            if looks_like_trna(p) and p.auth_asym_id not in wall:
                reasons[p.auth_asym_id] = "tRNA (auto-detected)"

    return {
        "wall_chain_ids": wall,
        "excluded_chain_ids": set(reasons.keys()),
        "reasons": reasons,
    }
