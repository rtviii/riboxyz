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
