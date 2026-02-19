# ribctl/lib/npet2/adapters/riboxyz_providers.py
"""
Providers that bridge riboxyz internals -> npet2 interfaces.

These are only usable inside the riboxyz repo where RibosomeOps, AssetType,
PTC_location, get_constriction are available. For standalone use, see
standalone_providers.py.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from npet2.core.ribosome_types import (
    RibosomeProfile,
    ProteinEntry,
    RNAEntry,
    PolymerEntry,
    AssemblyInfo,
    AssemblyPolymerInstance,
    AssemblyNonpolymerInstance,
)


def _convert_profile(ribo_struct) -> RibosomeProfile:
    """Convert a riboxyz RibosomeStructure to the npet2-internal RibosomeProfile."""

    def _nomenclature_strings(poly) -> list[str]:
        return [n.value if hasattr(n, "value") else str(n) for n in (poly.nomenclature or [])]

    proteins = [
        ProteinEntry(
            auth_asym_id=p.auth_asym_id,
            assembly_id=p.assembly_id,
            nomenclature=_nomenclature_strings(p),
            rcsb_pdbx_description=p.rcsb_pdbx_description,
            entity_poly_seq_length=p.entity_poly_seq_length,
        )
        for p in (ribo_struct.proteins or [])
    ]
    rnas = [
        RNAEntry(
            auth_asym_id=r.auth_asym_id,
            assembly_id=r.assembly_id,
            nomenclature=_nomenclature_strings(r),
            rcsb_pdbx_description=r.rcsb_pdbx_description,
            entity_poly_seq_length=r.entity_poly_seq_length,
        )
        for r in (ribo_struct.rnas or [])
    ]
    others = [
        PolymerEntry(
            auth_asym_id=o.auth_asym_id,
            assembly_id=o.assembly_id,
            nomenclature=_nomenclature_strings(o),
            rcsb_pdbx_description=o.rcsb_pdbx_description,
            entity_poly_seq_length=o.entity_poly_seq_length,
        )
        for o in (ribo_struct.other_polymers or [])
    ]

    assembly_map = None
    if ribo_struct.assembly_map:
        assembly_map = []
        for asm in ribo_struct.assembly_map:
            polys = [
                AssemblyPolymerInstance(
                    entity_id=inst.rcsb_polymer_entity_instance_container_identifiers.entity_id,
                    auth_asym_id=inst.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id,
                )
                for inst in (asm.polymer_entity_instances or [])
            ]
            nonpolys = None
            if asm.nonpolymer_entity_instances:
                nonpolys = [
                    AssemblyNonpolymerInstance(
                        entity_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.entity_id,
                        auth_asym_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.auth_asym_id,
                        auth_seq_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.auth_seq_id,
                    )
                    for inst in asm.nonpolymer_entity_instances
                ]
            assembly_map.append(AssemblyInfo(
                rcsb_id=asm.rcsb_id,
                polymer_entity_instances=polys,
                nonpolymer_entity_instances=nonpolys,
            ))

    return RibosomeProfile(
        rcsb_id=ribo_struct.rcsb_id,
        mitochondrial=bool(getattr(ribo_struct, "mitochondrial", False)),
        proteins=proteins,
        rnas=rnas,
        other_polymers=others,
        assembly_map=assembly_map,
    )


class RiboxyzStructureProvider:
    """Loads atoms + profile from local riboxyz assets (RibosomeOps)."""

    def fingerprint(self, rcsb_id: str) -> str:
        from ribctl.asset_manager.asset_types import AssetType
        p = AssetType.MMCIF.get_path(rcsb_id)
        return f"mmcif:{p}"

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        from ribctl.ribosome_ops import RibosomeOps
        from ribctl.asset_manager.asset_types import AssetType

        ro = RibosomeOps(rcsb_id)
        structure = ro.assets.biopython_structure()
        atoms = list(structure[0].get_atoms())
        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray([getattr(a, "element", "") or a.get_id()[0] for a in atoms])

        profile = _convert_profile(ro.profile)

        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(AssetType.MMCIF.get_path(rcsb_id)),
            "profile": profile,
            # Keep ro around for legacy stages that need biopython_structure etc.
            "ro": ro,
        }


class RiboxyzLandmarkProvider:
    def fingerprint(self, rcsb_id: str) -> str:
        return "ptc_via_trna+constriction_site:v1"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        from ribctl.lib.landmarks.ptc_via_trna import PTC_location
        from ribctl.lib.landmarks.constriction_site import get_constriction

        ptc = np.array(PTC_location(rcsb_id).location, dtype=np.float32)
        constr = np.array(get_constriction(rcsb_id), dtype=np.float32)
        return {"ptc_xyz": ptc, "constriction_xyz": constr}
