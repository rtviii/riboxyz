from typing import NewType, TypedDict
import typing
from ninja import Schema
from pydantic import BaseModel, create_model
from ribctl.lib.ribosome_types.types_ribosome import NonpolymericLigand, PolynucleotideClass, Protein, CytosolicProteinClass, RibosomeStructure
"""This file documents the possible requests that the API can receive."""


class ExogenousRNAByStruct(Schema):
    struct: str
    rnas: list[str]

class LigandByStructInstance(Schema):
     chemid: str
     name:str
     number:int

class LigandsByStruct(Schema): 
      title                  : str
      struct                 : str
      organism               : list[str]
      taxid                  : list[int]
      ligands                : list[LigandByStructInstance]


class PresentInStruct(Schema): 
      citation_title         : str
      description            : str
      expMethod              : str
      rcsb_id                : str
      resolution             : float
      src_organism_ids       : list[int]
    
class LigandInstance(Schema):
    # TODO: Huge todo. Refactor ligands into classes at the database level.
    chemicalId : str
    description: str
    polymer    : typing.Literal[False]
    presentIn  : PresentInStruct

class LigandlikeInstance(Schema):
    # TODO: See if ligandlike can be made into classes at the database level.
    polymer    : typing.Literal[True]
    description: str
    presentIn  : PresentInStruct

class PolymerMinimal(Schema):
    nomenclature: list[str]
    auth_asym_id: str
    entity_poly_seq_one_letter_code: str


class RibosomeHeader(Schema):
    rcsb_id:    str
    expMethod:  str
    resolution: float

    pdbx_keywords:      str | None
    pdbx_keywords_text: str | None

    rcsb_external_ref_id: list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year        : None | int
    citation_rcsb_authors: None | list[str]
    citation_title       : None | str
    citation_pdbx_doi    : None | str

    src_organism_ids: list[int]
    src_organism_names: list[str]

    host_organism_ids: list[int]
    host_organism_names: list[str]

class NeoStruct(Schema): 
      struct           : RibosomeHeader
      ligands          : list[str] | None
      rps              : list[PolymerMinimal]
      rnas             : list[PolymerMinimal] | None



class Residue(Schema):
    residue_name: str
    residue_id: int
    parent_auth_asym_id: str


class ExogenousRNAInStruct(Schema):
    struct: str
    rnas: list[str]

class NomenclatureClassMember(Schema):
    #TODO: Deprecate. lacks  ligand_like and host organism info. no surface ratio exists anymore. lazy design
    parent_resolution                  : float
    parent_year                        : int | None
    parent_method                      : str
    parent_citation                    :str
    parent_rcsb_id   : str | None

    pfam_accessions  : list[str]  | None
    pfam_comments    : list[str]  | None
    pfam_descriptions: list[str]  | None

    asym_ids:list[str]
    auth_asym_id: str

    src_organism_names : list[str]
    src_organism_ids   : list[int]

    uniprot_accession: list[str]| None
    rcsb_pdbx_description: str | None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str
    
    nomenclature:  list
    ligand_like:   bool | None



class BindingSiteChain(Schema): 
      sequence                : str
      nomenclature            : list[str]
      asym_ids                : list[str]
      residues                : list[Residue]

# export type LigandBindingSite = {
#   [ chainname    :string ] :
#   { sequence     : string
#     nomenclature : string[ ]
#     asym_ids     : string[ ]
#     residues     : Residue[]
#   }
# }

# ? https://stackoverflow.com/questions/72268685/pydantic-checks-on-newtype
class LigandBindingSite(Schema):
    __root__: dict[str, BindingSiteChain]
    def __getattr__(self, attr):
        return self.__root__[attr]



class PredictionSource(Schema):
    src: str
    auth_asym_id: str
    src_ids: list[int]
    
class PredictionTarget(Schema):
    tgt: str
    auth_asym_id: str
    tgt_ids: list[int]

class Alignement(Schema):
    src_aln: str
    tgt_aln: str
    aln_ids: list[int]

# export type LigandPrediction = {
#   [ polypeptide_class :string ] :
#   {
#   source   : {src    :string, auth_asym_id :string, src_ids:number[]},
#   target   : {tgt    :string, auth_asym_id :string, tgt_ids:number[]},
#   alignment: {src_aln:string, tgt_aln:string, aln_ids:number[]},
#   }
# }


class MixedLigand(Schema): 
      category           : str | None
      polymer            : bool
      description        : str
      chemicalId         : str |None
      present_in         : LigandBindingSite

# export interface MixedLigand{
#     category     ?: string,
#     polymer       : boolean,
#     description   : string;
#     chemicalId   ?: string
#     present_in: BindingSite
# }

class LigandClass(Schema):
    __root__: dict[str, list[MixedLigand]]
    def __getattr__(self, attr):
        return self.__root__[attr]


# export type LigandClass = {
#   [ligand_description:string]: MixedLigand[]
# }




# export type StructureBindingSites  =  {
#   [rcsb_id:string]: BindingSite[]
# }

class BindingSite(Schema):
    src_organism_ids: list[int]
    description: str
    citation_title: str
    auth_asym_id: str | None
    expMethod: str
    rcsb_id: str
    resolution: float

class StructureBindingSites(Schema):
    __root__: dict[str, list[BindingSite]]
    def __getattr__(self, attr):
        return self.__root__[attr]

# export type BindingSite  =  {
#                                  src_organism_ids   : number[],
#                                  description        : string,
#                                  citation_title     : string,
#                                  auth_asym_id     ? : string;    // if it's a polymer, i suppose.
#                                  expMethod          : string,
#                                  rcsb_id            : string,
#                                  resolution         : number,
# }
class NeoHomolog(Schema):
    parent: str
    orgname: list[str]
    orgid: list[int]
    protein: Protein
    title: str
# export interface NeoHomolog {
#   parent : string;
#   orgname: string[]
#   orgid  : number[]
#   protein: Protein;
#   title  : string
# }



class ProteinProfile(Schema):

    asym_ids: list[str]
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

    nomenclature:  list[CytosolicProteinClass]

    parent_resolution: float
    parent_year      : int
    parent_method    : str




class RNA(Schema):

    asym_ids: list[str]

    auth_asym_id: str
    nomenclature: list[PolynucleotideClass]
    parent_rcsb_id: str

    src_organism_names: list[str]
    host_organism_names: list[str]
    src_organism_ids: list[int]
    host_organism_ids: list[int]

    rcsb_pdbx_description: str | None

    entity_poly_strand_id: str
    entity_poly_seq_one_letter_code: str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length: int
    entity_poly_polymer_type: str
    entity_poly_entity_type: str

    ligand_like: bool

    parent_year      : int
    parent_resolution: float
    parent_method    : str



class StructureWithLigand(Schema):

    src_organism_ids: list[int]
    src_organism_names: list[str]
    rcsb_id: str
    expMethod: str
    resolution: float
    cryoem_exp_resolution: float | None
    citation_title: str
    pdbx_keywords_text: str | None


class LigandResponseShape(Schema):
    ligand: NonpolymericLigand
    presentIn: list[StructureWithLigand]


class BanClassMetadata(Schema):
    banClass: CytosolicProteinClass
    organisms: list[int]
    comments: list[list[str]]
    structs: list[str]

class RPSummary(Schema):
    organism_desc: list[str]
    organism_id  : list[int]
    uniprot      : list[str]
    parent       : str
    parent_reso  : float
    strand_id    : str
                        
# #TODO: Merge or disjoin with the BanClassMetadata class
class NomenclatureClass(Schema): 
      structs                  : list[str]
      rps                      : list[RPSummary]
      banClass                 : CytosolicProteinClass