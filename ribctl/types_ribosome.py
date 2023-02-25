import typing
from pydantic import BaseModel
from pydantic.dataclasses import dataclass
import dataclasses
import json






SSU_Proteins = typing.Literal["bS1" ,
"eS1" ,
"uS2" ,
"uS3" ,
"uS4" ,
"eS4" ,
"uS5" ,
"bS6" ,
"eS6" ,
"uS7" ,
"eS7" ,
"uS8" ,
"eS8" ,
"uS9" ,
"uS10",
"eS10",
"uS11",
"uS12",
"eS12",
"uS13",
"uS14",
"uS15",
"bS16",
"uS17",
"eS17",
"bS18",
"uS19",
"eS19",
"bS20",
"bS21",
"bTHX",
"eS21",
"eS24",
"eS25",
"eS26",
"eS27",
"eS28",
"eS30",
"eS31",
"RACK1"]



LSU_Proteins = typing.Literal["uL1",
"uL2",
"uL3",
"uL4",
"uL5",
"uL6",
"eL6",
"eL8",
"bL9",
"uL10",
"uL11",
"bL12",
"uL13",
"eL13",
"uL14",
"eL14",
"uL15",
"eL15",
"uL16",
"bL17",
"uL18",
"eL18",
"bL19",
"eL19",
"bL20",
"eL20",
"bL21",
"eL21",
"uL22",
"eL22",
"uL23",
"uL24",
"eL24",
"bL25",
"bL27",
"eL27",
"bL28",
"eL28",
"uL29",
"eL29",
"uL30",
"eL30",
"bL31",
"eL31",
"bL32",
"eL32",
"bL33",
"eL33",
"bL34",
"eL34",
"bL35",
"bL36",
"eL36",
"eL37",
"eL38",
"eL39",
"eL40",
"eL41",
"eL42",
"eL43",
"P1/P2" ]











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

