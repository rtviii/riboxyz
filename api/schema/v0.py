from typing import NewType, TypedDict
from ninja import Schema
from ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import Ligand, Protein, ProteinClass, RibosomeStructure
"""This file documents the possible requests that the API can receive."""



class PolymerMinimal(Schema):
    nomenclature: list[str]
    auth_asym_id: str
    entity_poly_seq_one_letter_code: str

class NeoStruct(Schema): 
      struct           : RibosomeStructure
      ligands          : list[str]
      rps              : list[PolymerMinimal]
      rnas             : list[PolymerMinimal]


class Residue(Schema):
    residue_name: str
    residue_id: int
    parent_auth_asym_id: str




class BindingSiteChain(Schema):
    sequence:str
    nomenclature:list[str]
    asym_ids:list[str]
    residues:list[Residue]
# export type LigandBindingSite = {
#   [ chainname    :string ] :
#   { sequence     : string
#     nomenclature : string[ ]
#     asym_ids     : string[ ]
#     residues     : Residue[]
#   }
# }
# TODO: Does this work?
# ? https://stackoverflow.com/questions/72268685/pydantic-checks-on-newtype
class LigandBindingSite(Schema):
    __root__: dict[str, BindingSiteChain]




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
class LigandPrediction(Schema):
    __root__: dict[str,dict[str, LigandBindingSite]]
# export type LigandPrediction = {
#   [ polypeptide_class :string ] :
#   {
#   source   : {src    :string, auth_asym_id :string, src_ids:number[]},
#   target   : {tgt    :string, auth_asym_id :string, tgt_ids:number[]},
#   alignment: {src_aln:string, tgt_aln:string, aln_ids:number[]},
#   }
# }


class MixedLigand(Schema):
    category: str | None
    polymer: bool
    description: str
    chemicalId: str |None
    present_in: LigandBindingSite

# export interface MixedLigand{
#     category     ?: string,
#     polymer       : boolean,
#     description   : string;
#     chemicalId   ?: string
#     present_in: BindingSite
# }

class LigandClass(Schema):
    __root__: dict[str, list[MixedLigand]]


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

    nomenclature:  list[ProteinClass]

    parent_resolution: float
    parent_year      : int
    parent_method    : str




class RNA(Schema):

    asym_ids: list[str]

    auth_asym_id: str
    nomenclature: list[RNAClass]
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
    ligand: Ligand
    presentIn: list[StructureWithLigand]


class BanClassMetadata(Schema):
    banClass: ProteinClass
    organisms: list[int]
    comments: list[list[str]]
    structs: list[RibosomeStructure]