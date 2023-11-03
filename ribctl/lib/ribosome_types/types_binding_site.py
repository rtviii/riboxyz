import json
import os
import typing
import pydantic
from pydantic import BaseModel
from Bio.PDB.Residue import Residue
from ribctl.lib.ribosome_types.types_ribosome import Polymer, PolymerClass

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
            seqid=seqid,
            resname=r.get_resname(),
            parent_auth_asym_id=chain_id,
            full_id=r.get_full_id()
        )

class BindingSiteChain(Polymer):
    residues: list[ResidueSummary]

class BindingSite(BaseModel):
    __root__: typing.Dict[str, BindingSiteChain]

    @staticmethod
    def path_nonpoly_ligand( rcsb_id: str, class_: str):
        RIBETL_DATA = os.environ.get('RIBETL_DATA')
        return os.path.join( str(RIBETL_DATA), rcsb_id.upper() , "LIGAND_" + class_.replace(" ", "_").upper() + ".json")

    # @staticmethod
    # def path_poly_factor( rcsb_id: str, class_: str, auth_asym_id: str):
    #     RIBETL_DATA = os.environ.get('RIBETL_DATA')
    #     return os.path.join( str(RIBETL_DATA), rcsb_id.upper(), "POLYMER_"+auth_asym_id +  "_" + class_.replace(" ", "_").upper() + ".json" )

    def save(self, filename: str):
        with open(filename, 'w') as outfile:
            json.dump(json.loads(self.json()), outfile, indent=4)
            print("Saved: ",filename)

    def __getattr__(self, attr):
        return super().dict()['__root__'].__getattribute__(attr)

    def __getitem__(self, attr):
        return super().dict()['__root__'].__getitem__(attr)

    def dict(self,):
        return super().dict()['__root__']

class PredictedResiduesPolymer(BaseModel):

    class PredictionSource(BaseModel):
        src: str
        src_ids: list[int]
        auth_asym_id: str

    class PredictionTarget(BaseModel):
        tgt: str
        tgt_ids: list[int]
        auth_asym_id: str

    class PredictionAlignments(BaseModel):
        aln_ids: list[int]
        src_aln: str
        tgt_aln: str

    source: PredictionSource
    target: PredictionTarget
    alignment: PredictionAlignments

class LigandPrediction(BaseModel):
    __root__: typing.Dict[PolymerClass, PredictedResiduesPolymer]

    def __getattr__(self, attr):
        return super().dict()['__root__'].__getattribute__(attr)

    def __getitem__(self, attr):
        return super().dict()['__root__'].__getitem__(attr)

    def dict(self,):
        return super().dict()['__root__']

# So, optimal transport for types.
