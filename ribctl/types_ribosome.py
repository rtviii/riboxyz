import typing
from pydantic import BaseModel
from pydantic.dataclasses import dataclass
import dataclasses
import json
from types_proteins import LSU_Proteins, SSU_Proteins




ProteinClass = LSU_Proteins | SSU_Proteins
RNAClass            = typing.Literal["5SrRNA" , "5.8SrRNA" , "12SrRNA" , "16SrRNA" , "21SrRNA" , "23SrRNA" , "25SrRNA" , "28SrRNA" , "35SrRNA" , "mRNA" , "tRNA"]
ProteinClassesArray = typing.get_args(ProteinClass)

class LastUpdate(BaseModel):
    date: str
    added_structure: str


class Protein(BaseModel):
    asym_ids    : list[str]
    auth_asym_id: str

    parent_rcsb_id   : str
    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]

    src_organism_names : list[str]
    host_organism_names: list[str]
    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    ligand_like: bool

    uniprot_accession: list[str]

    rcsb_pdbx_description: str | None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    nomenclature:  list[ProteinClass]



class RNA(BaseModel):
    asym_ids: list[str]

    auth_asym_id  : str
    nomenclature  : list[RNAClass]
    parent_rcsb_id: str

    src_organism_names : list[str]
    host_organism_names: list[str]
    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description: str | None
    # entity_polymer
    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    ligand_like: bool


class Ligand(BaseModel):
    chemicalId: str
    chemicalName: str
    formula_weight: float
    pdbx_description: str
    number_of_instances: int


class RibosomeStructure(BaseModel):
    rcsb_id: str
    expMethod: str
    resolution: float

    pdbx_keywords     : str | None
    pdbx_keywords_text: str | None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    # citation
    citation_year: int
    citation_rcsb_authors: list[str]
    citation_title: str
    citation_pdbx_doi: str
    # keywords
    # custom
    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    proteins: list[Protein]
    rnas    : list[RNA] | None
    ligands : list[Ligand] | None


#※--------------------------------------------------------

class InterProFamily(BaseModel):
    family_id: str
    type: str
    description: str

class GOClass(BaseModel):
    class_id: str
    annotation: str

class PFAMFamily(BaseModel): 
      family_id            : str
      annotation           : str
      family_type          : str

class NomeclatureClass(BaseModel):
    class_id:  ProteinClass | RNAClass

