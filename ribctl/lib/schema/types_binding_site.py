import json
import os
import typing
import pydantic
from pydantic import BaseModel, RootModel, field_serializer
from Bio.PDB.Residue import Residue
from ribctl.lib.schema.types_ribosome import Polymer  ,ResidueSummary
from ribctl.lib.types.polymer import PolymerClass

class BindingSiteChain(Polymer):
    bound_residues: list[ResidueSummary]

class BindingSite(BaseModel):
    source: str
    ligand: str
    radius: float
    chains: list[BindingSiteChain]

#TODO: The following four classes can be collapsed into one
class PredictionSource(BaseModel):
    source_seq    : str
    source_bound_residues: list[ResidueSummary]
    auth_asym_id  : str

#TODO: these two should just b ea single class tagged extending Polymer with "TARGET"| "SOURCE" and "bound residues"
class PredictionTarget(BaseModel):
    target_seq           : str
    target_bound_residues: list[ResidueSummary]
    auth_asym_id         : str

class PredictionAlignments(BaseModel):
    aligned_ids       : list[int]
    source_seq_aligned: str
    target_seq_aligned: str

class ResiduesMapping(BaseModel):

    polymer_class: PolymerClass # type: ignore
    source       : PredictionSource
    target       : PredictionTarget


class LigandTransposition(BaseModel):

    source                : str
    target                : str
    constituent_chains    : list[ResiduesMapping]
    purported_binding_site: BindingSite

    def save(self, filename: str):
        with open(filename, "w") as outfile:
            json.dump(self.model_dump(), outfile, indent=4)
            print("Saved: ", filename)
