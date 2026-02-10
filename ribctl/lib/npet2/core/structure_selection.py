# ribctl/lib/npet2/core/structure_selection.py
from __future__ import annotations

from typing import Iterable, Set

_TRNA_HINTS = ("trna", "transfer rna")


def _poly_description(poly) -> str:
    return (getattr(poly, "rcsb_pdbx_description", None) or "").strip()


def _poly_nomenclature_names(poly) -> list[str]:
    noms = getattr(poly, "nomenclature", None) or []
    out: list[str] = []
    for x in noms:
        # PolymerClass enum objects often have .name
        out.append(getattr(x, "name", str(x)))
    return out


def looks_like_trna(poly) -> bool:
    """
    Heuristic: profile-provided description/nomenclature.
    (You said you have a better classifier — you can replace this logic later.)
    """
    desc = _poly_description(poly).lower()
    if any(h in desc for h in _TRNA_HINTS):
        return True

    noms = " ".join(n.lower() for n in _poly_nomenclature_names(poly))
    if "trna" in noms:
        return True

    return False


def ribosome_wall_auth_asym_ids(
    profile,
    *,
    exclude_trna: bool = True,
    extra_exclude: Iterable[str] = (),
) -> Set[str]:
    """
    Define the tunnel 'wall' chains.

    Core rule:
      wall = proteins + rRNAs from the ribosome profile

    Benefits:
      - excludes HOH/ions/ligands/nonpolymers automatically (different chains/entities)
      - includes modified residues inside ribosomal polymers automatically
    """
    proteins = getattr(profile, "proteins", None) or []
    rnas = getattr(profile, "rnas", None) or []

    wall = {p.auth_asym_id for p in proteins} | {r.auth_asym_id for r in rnas}

    if exclude_trna:
        # Defensive: in some datasets a tRNA might be misfiled under rnas/other_polymers.
        others = getattr(profile, "other_polymers", None) or []
        all_polys = list(proteins) + list(rnas) + list(others)
        trna_ids = {p.auth_asym_id for p in all_polys if looks_like_trna(p)}
        wall -= trna_ids

    wall -= set(extra_exclude)
    return wall


def intersect_with_first_assembly(ro, chain_ids: Set[str]) -> Set[str]:
    """
    Optional safety: only use chains present in the first assembly.
    """
    try:
        asm = set(ro.first_assembly_auth_asym_ids())
        return chain_ids & asm
    except Exception:
        return chain_ids

# ribctl/lib/npet2/core/structure_selection.py  -- add this function:

def atom_inclusion_policy(profile, config, rcsb_id: str, ro) -> dict:
    """
    Central policy for which atoms go into occupancy calculations.

    INCLUDED (define tunnel walls):
      - Ribosomal proteins (all atoms, including modified residues)
      - rRNAs (all atoms, including modified nucleotides)

    EXCLUDED (treated as void / not wall):
      - Water molecules (HOH) -- solvent, should be inside tunnel
      - Ions (Mg2+, K+, etc.) -- solvent-associated
      - Nonpolymer ligands (antibiotics, spermidine, paromomycin, etc.)
      - tRNAs (configurable, default: excluded because they block PTC region)
      - Manually specified chains (config.occupancy_exclude_auth_asym_ids)
      - Known tunnel-debris chains (hardcoded per structure)

    Note: modified residues within ribosomal polymers ARE included because
    they are covalently part of the wall (e.g., pseudouridine, methylated bases).
    Waters/ions/ligands are on separate mmCIF entities/chains and are excluded
    by virtue of only selecting protein + rRNA auth_asym_ids.

    Returns dict with:
      - wall_chain_ids: set of auth_asym_ids for occupancy
      - excluded_chain_ids: set of auth_asym_ids that were explicitly removed
      - reason: dict mapping excluded chain_id -> reason string
    """
    from ribctl.lib.npet2.stages.legacy_minimal import _tunnel_debris_chains

    debris = _tunnel_debris_chains(rcsb_id, ro, profile)
    manual_exclude = list(getattr(config, "occupancy_exclude_auth_asym_ids", ()))
    exclude_trna = bool(getattr(config, "occupancy_exclude_trna", True))

    all_exclude = list(dict.fromkeys(debris + manual_exclude))

    wall = ribosome_wall_auth_asym_ids(
        profile,
        exclude_trna=exclude_trna,
        extra_exclude=all_exclude,
    )
    wall = intersect_with_first_assembly(ro, wall)

    # Build reason map for logging
    reasons = {}
    for c in debris:
        reasons[c] = "tunnel_debris (hardcoded)"
    for c in manual_exclude:
        if c not in reasons:
            reasons[c] = "config exclude"

    if exclude_trna:
        others = getattr(profile, "other_polymers", None) or []
        proteins = getattr(profile, "proteins", None) or []
        rnas = getattr(profile, "rnas", None) or []
        all_polys = list(proteins) + list(rnas) + list(others)
        for p in all_polys:
            if looks_like_trna(p) and p.auth_asym_id not in wall:
                reasons[p.auth_asym_id] = "tRNA (auto-detected)"

    return {
        "wall_chain_ids": wall,
        "excluded_chain_ids": set(reasons.keys()),
        "reasons": reasons,
    }