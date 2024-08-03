import json
import os
import typing
import pydantic
from pydantic import BaseModel, RootModel
from Bio.PDB.Residue import Residue
from ribctl.lib.schema.types_ribosome import Polymer, PolymerClass, PolynucleotideClass

AMINO_ACIDS = {
    "ALA": {"one_letter_code": "A", "charge": 0},
    "ARG": {"one_letter_code": "R", "charge": 1},
    "ASN": {"one_letter_code": "N", "charge": 0},
    "ASP": {"one_letter_code": "D", "charge": -1},
    "CYS": {"one_letter_code": "C", "charge": 0},
    "GLU": {"one_letter_code": "E", "charge": -1},
    "GLN": {"one_letter_code": "Q", "charge": 0},
    "GLY": {"one_letter_code": "G", "charge": 0},
    "HIS": {"one_letter_code": "H", "charge": 0},
    "ILE": {"one_letter_code": "I", "charge": 0},
    "LEU": {"one_letter_code": "L", "charge": 0},
    "LYS": {"one_letter_code": "K", "charge": 1},
    "MET": {"one_letter_code": "M", "charge": 0},
    "PHE": {"one_letter_code": "F", "charge": 0},
    "PRO": {"one_letter_code": "P", "charge": 0},
    "SER": {"one_letter_code": "S", "charge": 0},
    "THR": {"one_letter_code": "T", "charge": 0},
    "TRP": {"one_letter_code": "W", "charge": 0},
    "TYR": {"one_letter_code": "Y", "charge": 0},
    "VAL": {"one_letter_code": "V", "charge": 0},
}
NUCLEOTIDES = ["A", "T", "C", "G", "U"]


class ResidueSummary(BaseModel): 

    full_id            : typing.Optional[tuple[str, int, str, tuple[str, int, str]]]
    resname            : str
    seqid              : int
    parent_auth_asym_id: str

    def __hash__(self):
        return hash(
            self.get_resname() + str(self.get_seqid()) + self.get_parent_auth_asym_id()
        )

    def get_resname(self):
        return self.resname

    def get_seqid(self):
        (structure_id, model_id, chain_id, _) = self.full_id
        (hetero, seqid, insertion_code) = _
        return seqid

    def get_parent_auth_asym_id(self):
        (structure_id, assembly_id, chain_id, _) = self.full_id
        return chain_id

    @staticmethod
    def from_biopython_residue(r: Residue):

        (structure_id, model_id, chain_id, _) = r.get_full_id()
        (hetero, seqid, insertion_code) = _

        return ResidueSummary(
            seqid=seqid,
            resname=r.get_resname(),
            parent_auth_asym_id=chain_id,
            full_id=r.get_full_id(),
        )


class BindingSiteChain(Polymer):
    residues: list[ResidueSummary]


type BindingSite = list[BindingSiteChain]


class PredictionSource(BaseModel):

    source_seq: str
    source_seq_ids: list[int]
    auth_asym_id: str


class PredictionTarget(BaseModel):

    target_seq: str
    target_seq_ids: list[int]
    auth_asym_id: str


class PredictionAlignments(BaseModel):
    aligned_ids: list[int]
    source_seq_aligned: str
    target_seq_aligned: str


class PredictedResiduesPolymer(BaseModel):

    polymer_class: PolymerClass
    source: PredictionSource
    target: PredictionTarget
    alignment: PredictionAlignments


class LigandTransposition(BaseModel):
    source: str
    target: str
    constituent_chains: list[PredictedResiduesPolymer]
    purported_binding_site: BindingSite

    def save(self, filename: str):
        with open(filename, "w") as outfile:
            json.dump(self.model_dump(), outfile, indent=4)
            print("Saved: ", filename)
