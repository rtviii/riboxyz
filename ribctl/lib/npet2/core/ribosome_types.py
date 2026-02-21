# ribctl/lib/npet2/core/ribosome_types.py
"""
Standalone pydantic models for ribosome structure profiles as consumed by npet2.

These mirror the relevant subset of the riboxyz schema. When running inside the
riboxyz repo, the adapter converts full RibosomeStructure objects to these.
When running standalone, users provide profile JSON that validates against these
models directly.
"""
from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, field_validator

from .polymer_enum import PolymerClass, parse_polymer_class


# ---------------------------------------------------------------------------
# Polymer entries (slim -- only what the pipeline touches)
# ---------------------------------------------------------------------------

class PolymerEntry(BaseModel):
    """A single polymer chain as seen by the npet2 pipeline."""
    auth_asym_id          : str
    assembly_id           : int = 0
    nomenclature          : list[PolymerClass] = []
    rcsb_pdbx_description : Optional[str] = None
    entity_poly_seq_length: int = 0

    @field_validator("nomenclature", mode="before")
    @classmethod
    def _coerce_nomenclature(cls, v):
        """Accept raw strings (from JSON) and convert to enum members."""
        if not v:
            return []
        out = []
        for item in v:
            if isinstance(item, str):
                out.append(parse_polymer_class(item))
            else:
                out.append(item)
        return out


class ProteinEntry(PolymerEntry):
    """Protein chain. Inherits all pipeline-relevant fields from PolymerEntry."""
    pass


class RNAEntry(PolymerEntry):
    """RNA chain."""
    pass


# ---------------------------------------------------------------------------
# Assembly map (for first-assembly filtering)
# ---------------------------------------------------------------------------

class AssemblyPolymerInstance(BaseModel):
    entity_id: str
    auth_asym_id: str


class AssemblyNonpolymerInstance(BaseModel):
    entity_id: str
    auth_asym_id: str
    auth_seq_id: str = ""


class AssemblyInfo(BaseModel):
    rcsb_id: str
    polymer_entity_instances: list[AssemblyPolymerInstance] = []
    nonpolymer_entity_instances: Optional[list[AssemblyNonpolymerInstance]] = None


# ---------------------------------------------------------------------------
# Top-level profile
# ---------------------------------------------------------------------------

class RibosomeProfile(BaseModel):
    """
    Minimal ribosome structure profile for the npet2 pipeline.

    This is the only "ribosome knowledge" npet2 needs beyond the mmCIF itself
    and the landmark coordinates. It tells the pipeline which chains form the
    tunnel walls (proteins + rRNAs), which to exclude (tRNAs, factors), and
    whether the structure is mitochondrial.

    Can be loaded from:
      - riboxyz API: GET /structures/{rcsb_id}/profile
      - Local JSON file (--profile path/to/profile.json)
      - Constructed programmatically by the riboxyz adapter
    """
    rcsb_id: str
    mitochondrial: bool = False

    proteins: list[ProteinEntry] = []
    rnas: list[RNAEntry] = []
    other_polymers: list[PolymerEntry] = []

    assembly_map: Optional[list[AssemblyInfo]] = None

    def all_polymers(self) -> list[PolymerEntry]:
        return [*self.proteins, *self.rnas, *self.other_polymers]

    def first_assembly_auth_asym_ids(self) -> list[str]:
        """Return auth_asym_ids from the first assembly, or raise."""
        if not self.assembly_map:
            raise ValueError(
                f"No assembly_map in profile for {self.rcsb_id}. "
                "Provide one, or set occupancy_chain_mode='walls_only' to skip assembly filtering."
            )
        first = self.assembly_map[0]
        ids = [inst.auth_asym_id for inst in first.polymer_entity_instances]
        if first.nonpolymer_entity_instances:
            ids.extend(inst.auth_asym_id for inst in first.nonpolymer_entity_instances)
        return ids


# ---------------------------------------------------------------------------
# Landmark types
# ---------------------------------------------------------------------------

class PTCInfo(BaseModel):
    """Peptidyl transferase center location."""
    location: list[float]

    @field_validator("location")
    @classmethod
    def _validate_location(cls, v):
        if len(v) != 3:
            raise ValueError(f"PTC location must be [x, y, z], got {len(v)} values")
        return [float(x) for x in v]


class ConstrictionInfo(BaseModel):
    """Constriction site location."""
    location: list[float]

    @field_validator("location")
    @classmethod
    def _validate_location(cls, v):
        if len(v) != 3:
            raise ValueError(f"Constriction location must be [x, y, z], got {len(v)} values")
        return [float(x) for x in v]