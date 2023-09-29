from typing import Any, Optional
from enum import Enum, auto
import typing
import typing
from typing import NewType
from pydantic import BaseModel
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import LSUProteinClass, NonpolymericLigandClass, PolymericFactorClass, RNAClass, SSUProteinClass

RCSB_ID = NewType('RCSB_ID', str)

ProteinClass = typing.Union[LSUProteinClass , SSUProteinClass]
PolymerClass = typing.Union[ProteinClass, RNAClass]

class Polymer(BaseModel):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

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
        
     

        # nonpoly["nonpolymer_comp"]["chem_comp"]["id"]
        # nonpoly["nonpolymer_comp"]["chem_comp"]["name"]
        # nonpoly["nonpolymer_comp"]["chem_comp"]["three_letter_code"]

        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_container_identifiers"][
        #     "drugbank_id"
        # ]
        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_info"]["cas_number"]
        # nonpoly["nonpolymer_comp"]["drugbank"]["drugbank_info"]["description"]

        # nonpoly["nonpolymer_comp"]["rcsb_chem_comp_target"]

        # {
        #     "interaction_type": "target",
        #     "name": "Spermine synthase",
        #     "provenance_source": "DrugBank",
        #     "reference_database_accession_code": "P52788",
        #     "reference_database_name": "UniProt",
        # },

class PolymericFactor(Polymer): 
    def __hash__(self) -> int:
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    nomenclature: list[PolymericFactorClass] 
    
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

    pdbx_keywords:      Optional[str]
    pdbx_keywords_text: Optional[str]

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

    proteins            : list[Protein]
    rnas                : Optional[list[RNA] ]
    nonpolymeric_ligands: Optional[list[NonpolymericLigand] ]
    polymeric_factors   : Optional[list[PolymericFactor] ]
    # This includes DNA-RNA hybrid strands, DNA and all other polymers
    other_polymers   : Optional[list[Polymer]]
    
    @staticmethod
    def from_json_profile(d: Any):
        return RibosomeStructure(**d)

class ProteinClassEnum(Enum):
    bS1   = "bS1"
    eS1   = "eS1"
    uS2   = "uS2"
    uS3   = "uS3"
    uS4   = "uS4"
    eS4   = "eS4"
    uS5   = "uS5"
    bS6   = "bS6"
    eS6   = "eS6"
    uS7   = "uS7"
    eS7   = "eS7"
    uS8   = "uS8"
    eS8   = "eS8"
    uS9   = "uS9"
    uS10  = "uS10"
    eS10  = "eS10"
    uS11  = "uS11"
    uS12  = "uS12"
    eS12  = "eS12"
    uS13  = "uS13"
    uS14  = "uS14"
    uS15  = "uS15"
    bS16  = "bS16"
    uS17  = "uS17"
    eS17  = "eS17"
    bS18  = "bS18"
    uS19  = "uS19"
    eS19  = "eS19"
    bS20  = "bS20"
    bS21  = "bS21"
    bTHX  = "bTHX"
    eS21  = "eS21"
    eS24  = "eS24"
    eS25  = "eS25"
    eS26  = "eS26"
    eS27  = "eS27"
    eS28  = "eS28"
    eS30  = "eS30"
    eS31  = "eS31"
    RACK1 = "RACK1"
    uL1   = "uL1"
    uL2   = "uL2"
    uL3   = "uL3"
    uL4   = "uL4"
    uL5   = "uL5"
    uL6   = "uL6"
    eL6   = "eL6"
    eL8   = "eL8"
    bL9   = "bL9"
    uL10  = "uL10"
    uL11  = "uL11"
    bL12  = "bL12"
    uL13  = "uL13"
    eL13  = "eL13"
    uL14  = "uL14"
    eL14  = "eL14"
    uL15  = "uL15"
    eL15  = "eL15"
    uL16  = "uL16"
    bL17  = "bL17"
    uL18  = "uL18"
    eL18  = "eL18"
    bL19  = "bL19"
    eL19  = "eL19"
    bL20  = "bL20"
    eL20  = "eL20"
    bL21  = "bL21"
    eL21  = "eL21"
    uL22  = "uL22"
    eL22  = "eL22"
    uL23  = "uL23"
    uL24  = "uL24"
    eL24  = "eL24"
    bL25  = "bL25"
    bL27  = "bL27"
    eL27  = "eL27"
    bL28  = "bL28"
    eL28  = "eL28"
    uL29  = "uL29"
    eL29  = "eL29"
    uL30  = "uL30"
    eL30  = "eL30"
    bL31  = "bL31"
    eL31  = "eL31"
    bL32  = "bL32"
    eL32  = "eL32"
    bL33  = "bL33"
    eL33  = "eL33"
    bL34  = "bL34"
    eL34  = "eL34"
    bL35  = "bL35"
    bL36  = "bL36"
    eL36  = "eL36"
    eL37  = "eL37"
    eL38  = "eL38"
    eL39  = "eL39"
    eL40  = "eL40"
    eL41  = "eL41"
    eL42  = "eL42"
    eL43  = "eL43"
    P1P2  = "P1P2"

class RNAClassEnum(Enum):
    #TODO: Assembly missing

    mt_rRNA_12S = "mt12SrRNA" # mitochondrial
    mt_rRNA_16S = "mt16SrRNA" # mitochondrial

    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA" #  c-bacterial or mitochondrial
    rRNA_23S  = "23SrRNA" # bacterial
    rRNA_25S  = "25SrRNA" # plants

    rRNA_5_8S = "5.8SrRNA" # eukaryotic
    rRNA_18S  = "18SrRNA" # eukaryotic
    rRNA_28S  = "28SrRNA" # eukaryotic

class ElongationFactor(Enum):
    # Eukaryotic
    eEF1A = "eEF1A"
    eEF1B = "eEF1B"
    EFsec = "EFsec"
    eEF2  = "eEF2"
    mtEF4 = "mtEF4"
    eIF5A = "eIF5A"
    eEF3  = "eEF3"
    # Bacterial
    EF_Tu = "EF-Tu "
    EF_Ts = "EF-Ts"
    SelB  = "SelB"
    EF_G  = "EF-G"
    EF4   = "EF4"
    EF_P  = "EF-P"
    Tet_O = "Tet(O)"
    Tet_M = "Tet(M)"
    RelA  = "RelA"
    BipA  = "BipA"
    # Archaeal
    aEF1A = "aEF1A"
    aEF2  = "aEF2"
    aIF5A = "aIF5A"

PolymerClass_ = typing.Union[RNAClassEnum, ProteinClassEnum, PolymericFactorClassEnum]