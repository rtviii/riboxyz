from enum import Enum
import json
from typing import List, Literal, Optional, Union
from pydantic import BaseModel, field_serializer

from ribctl.lib.schema.types_ribosome import CytosolicProteinClass, MitochondrialProteinClass, PolynucleotideClass

raw_chain = {
            "assembly_id": 1,
            "asym_ids": [
                "B",
                "CC"
            ],
            "auth_asym_id": "S0",
            "parent_rcsb_id": "5DAT",
            "src_organism_names": [
                "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"
            ],
            "host_organism_names": [],
            "src_organism_ids": [
                559292
            ],
            "host_organism_ids": [],
            "rcsb_pdbx_description": "40S ribosomal protein S0-A",
            "entity_poly_strand_id": "S0,s0",
            "entity_poly_seq_one_letter_code": "SLPATFDLTPEDAQLLLAANTHLGARNVQVHQEPYVFNARPDGVHVINVGKTWEKLVLAARIIAAIPNPEDVVAISSRTFGQRAVLKFAAHTGATPIAGRFTPGSFTNYITRSFKEPRLVIVTDPRSDAQAIKEASYVNIPVIALTDLDSPSEFVDVAIPCNNRGKHSIGLIWYLLAREVLRLRGALVDRTQPWSIMPDLYFYRDPEEVEQQVAEEATTEEAGEEEAKEEVTEEQAEATEWAEENADNVEW",
            "entity_poly_seq_one_letter_code_can": "SLPATFDLTPEDAQLLLAANTHLGARNVQVHQEPYVFNARPDGVHVINVGKTWEKLVLAARIIAAIPNPEDVVAISSRTFGQRAVLKFAAHTGATPIAGRFTPGSFTNYITRSFKEPRLVIVTDPRSDAQAIKEASYVNIPVIALTDLDSPSEFVDVAIPCNNRGKHSIGLIWYLLAREVLRLRGALVDRTQPWSIMPDLYFYRDPEEVEQQVAEEATTEEAGEEEAKEEVTEEQAEATEWAEENADNVEW",
            "entity_poly_seq_length": 251,
            "entity_poly_polymer_type": "Protein",
            "entity_poly_entity_type": "polypeptide(L)",
            "nomenclature": [
                "uS2"
            ],
            "pfam_accessions": [
                "PF00318"
            ],
            "pfam_comments": [
                "NULL"
            ],
            "pfam_descriptions": [
                "Ribosomal protein S2"
            ],
            "uniprot_accession": [
                "P32905"
            ]
        }

class Polymer(BaseModel):

    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    # def to_dict(self):
    #     """A hack for enum.union to work with pydantic BaseModel. Otherwise EnumUnion instances are represented as <MitochondrialProteinClass.mL64: 'mL64'> etc.(Correct is "mL64")"""
    #     return json.loads(self.json())


    # FIXME:
    # def to_SeqRecord(self)->SeqRecord:
    #     return SeqRecord(self.entity_poly_seq_one_letter_code_can, id=f"{self.parent_rcsb_id}.{self.auth_asym_id}", 
    #                      description="{}|{}".format(self.src_organism_ids[0],self.rcsb_pdbx_description))
    def __str__(self) -> str:
        return super().model_dump_json()

    assembly_id: int

    asym_ids    : list[str]
    auth_asym_id: str

    parent_rcsb_id   : str

    src_organism_names : list[str]
    host_organism_names: list[str]

    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description: Optional[str] 

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    @field_serializer('nomenclature')
    def serialize_dt(self, nomenclature:list[PolynucleotideClass], _info):
        return [x.name for x in nomenclature]
    nomenclature:  list[PolynucleotideClass]

class Protein(Polymer):


    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]
    uniprot_accession: list[str]

    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    @staticmethod
    def from_polymer(p: Polymer,
                      **kwargs): 

        if ( kwargs["pfams"] != None and len(kwargs["pfams"]) > 0 ):
            pfam_comments     = list( set( [ pfam["rcsb_pfam_comment"    ] for pfam in kwargs["pfams"] ] ) )
            pfam_descriptions = list( set( [ pfam["rcsb_pfam_description"] for pfam in kwargs["pfams"] ] ) )
            pfam_accessions   = list( set( [ pfam["rcsb_pfam_accession"  ] for pfam in kwargs["pfams"] ] ) )

        else:
            pfam_comments     = []
            pfam_descriptions = []
            pfam_accessions   = []

        
        return Protein(**{
              **p.model_dump(),
              "pfam_accessions"   : pfam_accessions,
              "pfam_comments"     : pfam_comments,
              "pfam_descriptions" : pfam_descriptions,
              "uniprot_accession" : [ entry["rcsb_id"] for entry in kwargs["uniprots"] ] if kwargs["uniprots"] != None and len(kwargs["uniprots"]) > 0 else []
        })


    def to_polymer(self)->Polymer:
        return Polymer(**self.model_dump())


class RNA(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)


class NonpolymericLigand(BaseModel)  : 

      class NonpolymerComp(BaseModel):
        class ChemComp(BaseModel):
            class ChemCompContainerIdentifiers(BaseModel):
                id: str
                name: str
                three_letter_code: str
            chem_comp_container_identifiers: ChemCompContainerIdentifiers

        class Drugbank(BaseModel):
            class DrugbankInfo(BaseModel):
                cas_number: Optional[str]
                description: Optional[str]
            class DrugbankContainerIdentifiers(BaseModel):
                drugbank_id: str
            drugbank_container_identifiers: DrugbankContainerIdentifiers
            drugbank_info: DrugbankInfo

        class RcsbChemCompTarget(BaseModel):

            interaction_type                 : Optional[str] 
            name                             : Optional[str] 
            provenance_source                : Optional[str] 
            reference_database_accession_code: Optional[str] 
            reference_database_name          : Optional[str] 

        
        chemp_comp           : Optional[ChemComp]
        drugbank             : Optional[Drugbank]
        rcsb_chem_comp_target: Optional[list[RcsbChemCompTarget]]

      chemicalId         : str
      chemicalName       : str
      formula_weight     : Optional[float]

      pdbx_description   : str
      number_of_instances: int

      nonpolymer_comp: NonpolymerComp
        

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

    class NonpolymerEntityInstance(BaseModel):

        class NonpolymerEntityInstanceContainerIdentifiers(BaseModel):

            entity_id   : str
            auth_asym_id: str
            auth_seq_id : str

        rcsb_nonpolymer_entity_instance_container_identifiers: NonpolymerEntityInstanceContainerIdentifiers

    class PolymerEntityInstance(BaseModel):
        class PolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str

        rcsb_polymer_entity_instance_container_identifiers: PolymerEntityInstanceContainerIdentifiers

    rcsb_id                    : str # ex. 5AFI-1
    nonpolymer_entity_instances: Optional[list[NonpolymerEntityInstance]]
    polymer_entity_instances   : list[PolymerEntityInstance]


class RibosomeStructure(BaseModel):

    rcsb_id   : str
    expMethod : str
    resolution: float


    pdbx_keywords         : Optional[str]
    pdbx_keywords_text    : Optional[str]

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year        : Optional[int]
    citation_rcsb_authors: Optional[list[str]]
    citation_title       : Optional[str]
    citation_pdbx_doi    : Optional[str]

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    assembly_map: list[AssemblyInstancesMap]
    mitochondrial: bool

    # proteins            : list[Any]
    # rnas                : list[Any]
    # nonpolymeric_ligands: list[Any]
    # polymeric_factors   : list[Any]

    proteins            : list[Protein]
    rnas                : list[RNA]

    # polymeric_factors   : list[LifecycleFactor]
    #? This includes DNA-RNA hybrid strands, DNA and all other polymers
    other_polymers   : list[Polymer]

    nonpolymeric_ligands: list[NonpolymericLigand]

    
    # TODO: Deprecate this
    # @staticmethod
    # def from_json_profile(d: Any):
    #     return RibosomeStructure(**d)

