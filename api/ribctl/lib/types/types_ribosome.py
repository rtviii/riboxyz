from typing import Any, Optional
import typing
from typing import NewType
from pydantic import BaseModel
from api.ribctl.lib.types.types_poly_nonpoly_ligand import LSU_Proteins, NonpolymericLigandClass, PolymericFactorClass, RNAClass, SSU_Proteins

RCSB_ID = NewType('RCSB_ID', str)


ProteinClass = typing.Union[LSU_Proteins , SSU_Proteins]
PolymerClass = typing.Union[ProteinClass, RNAClass]

class Polymer(BaseModel):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    assembly_id: int

    asym_ids: list[str]
    auth_asym_id: str

    parent_rcsb_id   : str

    src_organism_names : list[str]
    host_organism_names: list[str]
    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description: str | None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    nomenclature:  list[PolymerClass]

class Protein(Polymer):

    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]

    uniprot_accession: list[str]

    def to_poly(self)->Polymer:
        return Polymer(**self.dict())

class RNA(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    pass


class NonpolymericLigand(BaseModel)  : 

      chemicalId         : str
      chemicalName       : str
      formula_weight     : None | float
      pdbx_description   : str
      number_of_instances: int
    #   nomenclature       : list[NonpolymericLigandClass]

class PolymericFactor(Polymer): 
    def __hash__(self) -> int:
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    nomenclature: list[PolymericFactorClass] 
    
class NonpolymerEntityInstance(BaseModel):
    class NonpolymerEntityInstanceContainerIdentifiers(BaseModel):
        entity_id: str
        auth_asym_id: str
        auth_seq_id: str

    rcsb_nonpolymer_entity_instance_container_identifiers: NonpolymerEntityInstanceContainerIdentifiers

class PolymerEntityInstance(BaseModel):
    class PolymerEntityInstanceContainerIdentifiers(BaseModel):
        entity_id: str
        auth_asym_id: str
    rcsb_polymer_entity_instance_container_identifiers: PolymerEntityInstanceContainerIdentifiers


class AssemblyInstancesMap(BaseModel):
    """
    This basically specifies which assembly an instnace of a polymer or a nonpolymer belongs to. 
    Certain PDB structures come with more than a single physical model/assembly packaged in the file,
    hence every chain and many ligands might be present in 2 or more instances. 
    
    The RNA/Protein/Ligand suffices to characterizes all instances, yet, to resolve duplicate chains 
    precisely in space, this information is needed.
    
    assemblies{
    rcsb_id 
   	nonpolymer_entity_instances{
  
      rcsb_nonpolymer_entity_instance_container_identifiers{
        entity_id
        comp_id
        auth_asym_id
        rcsb_id
        auth_seq_id
      }
    }
    polymer_entity_instances{
       rcsb_polymer_entity_instance_container_identifiers {  
        entity_id
        asym_id
        auth_asym_id
        entry_id
        entity_id
      }
    }
  }
    """
    rcsb_id                    : str # 5AFI-1
    nonpolymer_entity_instances: Optional[list[NonpolymerEntityInstance]]
    polymer_entity_instances   : list[PolymerEntityInstance]

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

    assembly_map: list[AssemblyInstancesMap]

    proteins            : list[Protein]
    rnas                : list[RNA] | None
    nonpolymeric_ligands: list[NonpolymericLigand] | None
    polymeric_factors   : list[PolymericFactor] | None
    
    @staticmethod
    def from_json_profile(d: Any):
        return RibosomeStructure(**d)