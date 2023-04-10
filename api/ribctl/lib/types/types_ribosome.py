from typing import Any
import typing
from pydantic import BaseModel
from .types_polymer import LSU_Proteins, RNAClass, SSU_Proteins

ProteinClass = typing.Union[LSU_Proteins , SSU_Proteins]
PolymerClass = typing.Union[ProteinClass, RNAClass]

class LastUpdate(BaseModel):
    date: str
    added_structure: str

class Polymer(BaseModel):
    asym_ids: list[str]
    auth_asym_id: str

    parent_rcsb_id   : str

    src_organism_names : list[str]
    host_organism_names: list[str]
    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    ligand_like: bool

    rcsb_pdbx_description: str | None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    nomenclature:  list[PolymerClass]

class Protein(Polymer):

    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]

    uniprot_accession: list[str]

    def to_poly(self)->Polymer:
        return Polymer(**self.dict())

class RNA(Polymer):
    pass

class Ligand(BaseModel)  : 

      chemicalId         : str
      chemicalName       : str
      formula_weight     : None | float
      pdbx_description   : str
      number_of_instances: int

class RibosomeStructure(BaseModel):

    rcsb_id   : str
    expMethod : str
    resolution: float

    pdbx_keywords:      str | None
    pdbx_keywords_text: str | None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year        : None | int
    citation_rcsb_authors: None | list[str]
    citation_title       : None | str
    citation_pdbx_doi    : None | str

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    proteins: list[Protein]
    rnas    : list[RNA] | None
    ligands : list[Ligand] | None
    



    @staticmethod
    def from_json_profile(d: Any):
        return RibosomeStructure(**d)

class InterProFamily(BaseModel):
    family_id: str
    type: str
    description: str

class GOClass(BaseModel):
    class_id: str
    annotation: str

class PFAMFamily(BaseModel):
    family_id: str
    annotation: str
    family_type: str