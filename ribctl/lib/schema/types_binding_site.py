import json
import os
import typing
import pydantic
from pydantic import BaseModel, RootModel
from Bio.PDB.Residue import Residue
from ribctl.lib.schema.types_ribosome import Polymer, PolymerClass, PolynucleotideClass

AMINO_ACIDS = {
    "ALA": 0,
    'ARG': 1,
    'ASN': 0,
    'ASP': -1,
    'CYS': 0,
    'GLN': 0,
    'GLU': -1,
    'GLY': 0,
    'HIS': 0,
    'ILE': 0,
    'LEU': 0,
    'LYS': 1,
    'MET': 0,
    'PHE': 0,
    'PRO': 0,
    'SER': 0,
    'THR': 0,
    'TRP': 0,
    'TYR': 0,
    'VAL': 0,
    'SEC': 0,
    'PYL': 0
}

NUCLEOTIDES = ['A', 'T', 'C', 'G', 'U']

class ResidueSummary(BaseModel):
    full_id            : tuple[str, int, str, tuple[str, int, str]]
    resname            : str
    seqid              : int
    parent_auth_asym_id: str

    def __hash__(self):
        return hash(self.get_resname() + str(self.get_seqid()) + self.get_parent_auth_asym_id())

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
            seqid               = seqid,
            resname             = r.get_resname(),
            parent_auth_asym_id = chain_id,
            full_id             = r.get_full_id()
        )

class BindingSiteChain(Polymer):
    residues: list[ResidueSummary]

type BindingSite = list[BindingSiteChain]




class PredictionSource(BaseModel): 

  source_seq                 : str
  source_seq_ids             : list[int]
  auth_asym_id               : str

class PredictionTarget(BaseModel):

    target_seq         : str
    target_seq_ids     : list[int]
    auth_asym_id: str

class PredictionAlignments(BaseModel):
    aligned_ids: list[int]
    source_seq_aligned: str
    target_seq_aligned: str
class PredictedResiduesPolymer(BaseModel):


    polymer_class         : PolymerClass
    source                : PredictionSource
    target                : PredictionTarget
    alignment             : PredictionAlignments

class LigandTransposition(BaseModel):
    source: str
    target: str
    constituent_chains    : list[PredictedResiduesPolymer]
    purported_binding_site: BindingSite

    def save(self, filename: str):
        with open(filename, 'w') as outfile:
            json.dump(self.model_dump(), outfile, indent=4)
            print("Saved: ",filename)

