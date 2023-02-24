import typing
from pydantic import BaseModel

ProteinClass = typing.Literal["eL39", "eL38", "eL37", "eL36", "eL34", "eL33", "eL32", "eL31", "eL30", "eL6", "uL2", "uL1", "uL4", "uL3", "uL6", "uL5", "bL36", "eL29", "bL35", "eL28", "uL30", "eL27", "bL32", "bL31", "eL24", "bL34", "bL33", "eL22", "eL21", "eL20", "eL19", "bL25", "eL18", "uL22", "bL27", "eL15", "bL21", "eL14", "bL20", "eL13", "uL29", "bL9", "uL23", "uL24", "bL28", "uL10", "uL11", "bL12",
    "uL18", "eL43", "uL16", "eL42", "eL41", "uL14", "eL40", "uL15", "uL13", "bL17", "bL19", "eS12", "eS10", "bS20", "eS31", "eS30", "uS19", "uS17", "uS15", "uS13", "uS14", "RACK1", "eS19", "eS17", "bS21", "eS24", "eS1", "eS21", "bS6", "eS4", "eS7", "uS11", "eS6", "uS12", "eS8", "uS10", "bTHX", "uS3", "uS2", "bS18", "uS5", "uS4", "uS7", "uS9", "uS8", "bS16", "eS28", "eS27", "eS26", "eS25", "AMBIGUOUS"]
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

