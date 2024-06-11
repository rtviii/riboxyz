import json
from typing import Dict, Optional
from enum import Enum
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, field_serializer
from ribctl.lib.enumunion import enum_union

# TODO:
# |********************************************************************************************************|
# | https://docs.google.com/spreadsheets/d/1mapshbn1ArofPN-Omu8GG5QdcwlJ0ym0BlN252kkUBU/edit#gid=815712128 |
# |********************************************************************************************************|
class tRNA(str,Enum):
    tRNA = "tRNA"

class MitochondrialProteinClass(str,Enum):
    # mSSU
    bS1m  = "bS1m"
    uS2m  = "uS2m"
    uS3m  = "uS3m"
    uS4m  = "uS4m"
    uS5m  = "uS5m"
    bS6m  = "bS6m"
    uS7m  = "uS7m"
    uS8m  = "uS8m"
    uS9m  = "uS9m"
    uS10m = "uS10m"
    uS11m = "uS11m"
    uS12m = "uS12m"
    uS13m = "uS13m"
    uS14m = "uS14m"
    uS15m = "uS15m"
    bS16m = "bS16m"
    uS17m = "uS17m"
    bS18m = "bS18m"
    uS19m = "uS19m"
    bS21m = "bS21m"
    mS22  = "mS22"
    mS23  = "mS23"
    mS25  = "mS25"
    mS26  = "mS26"
    mS27  = "mS27"
    mS29  = "mS29"
    mS31  = "mS31"
    mS33  = "mS33"
    mS34  = "mS34"
    mS35  = "mS35"
    mS37  = "mS37"
    mS38  = "mS38"
    mS39  = "mS39"
    mS40  = "mS40"
    mS41  = "mS41"
    mS42  = "mS42"
    mS43  = "mS43"
    mS44  = "mS44"
    mS45  = "mS45"
    mS46  = "mS46"
    mS47  = "mS47"

    # mLSU
    uL1m  = "uL1m"
    uL2m  = "uL2m"
    uL3m  = "uL3m"
    uL4m  = "uL4m"
    uL5m  = "uL5m"
    uL6m  = "uL6m"
    bL9m  = "bL9m"
    uL10m = "uL10m"
    uL11m = "uL11m"
    bL12m = "bL12m"
    uL13m = "uL13m"
    uL14m = "uL14m"
    uL15m = "uL15m"
    uL16m = "uL16m"
    bL17m = "bL17m"
    uL18m = "uL18m"
    bL19m = "bL19m"
    bL20m = "bL20m"
    bL21m = "bL21m"
    uL22m = "uL22m"
    uL23m = "uL23m"
    uL24m = "uL24m"
    bL27m = "bL27m"
    bL28m = "bL28m"
    uL29m = "uL29m"
    uL30m = "uL30m"
    bL31m = "bL31m"
    bL32m = "bL32m"
    bL33m = "bL33m"
    bL34m = "bL34m"
    bL35m = "bL35m"
    bL36m = "bL36m"
    mL37  = "mL37"
    mL38  = "mL38"
    mL39  = "mL39"
    mL40  = "mL40"
    mL41  = "mL41"
    mL42  = "mL42"
    mL43  = "mL43"
    mL44  = "mL44"
    mL45  = "mL45"
    mL46  = "mL46"
    mL48  = "mL48"
    mL49  = "mL49"
    mL50  = "mL50"
    mL51  = "mL51"
    mL52  = "mL52"
    mL53  = "mL53"
    mL54  = "mL54"
    mL57  = "mL57"
    mL58  = "mL58"
    mL59  = "mL59"
    mL60  = "mL60"
    mL61  = "mL61"
    mL62  = "mL62"
    mL63  = "mL63"
    mL64  = "mL64"
    mL65  = "mL65"
    mL66  = "mL66"
    mL67  = "mL67"

class CytosolicProteinClass(str,Enum):
    # SSU
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
    # LSU
    uL1  = "uL1"
    uL2  = "uL2"
    uL3  = "uL3"
    uL4  = "uL4"
    uL5  = "uL5"
    uL6  = "uL6"
    eL6  = "eL6"
    eL8  = "eL8"
    bL9  = "bL9"
    uL10 = "uL10"
    uL11 = "uL11"
    bL12 = "bL12"
    uL13 = "uL13"
    eL13 = "eL13"
    uL14 = "uL14"
    eL14 = "eL14"
    uL15 = "uL15"
    eL15 = "eL15"
    uL16 = "uL16"
    bL17 = "bL17"
    uL18 = "uL18"
    eL18 = "eL18"
    bL19 = "bL19"
    eL19 = "eL19"
    bL20 = "bL20"
    eL20 = "eL20"
    bL21 = "bL21"
    eL21 = "eL21"
    uL22 = "uL22"
    eL22 = "eL22"
    uL23 = "uL23"
    uL24 = "uL24"
    eL24 = "eL24"
    bL25 = "bL25"
    bL27 = "bL27"
    eL27 = "eL27"
    bL28 = "bL28"
    eL28 = "eL28"
    uL29 = "uL29"
    eL29 = "eL29"
    uL30 = "uL30"
    eL30 = "eL30"
    bL31 = "bL31"
    eL31 = "eL31"
    bL32 = "bL32"
    eL32 = "eL32"
    bL33 = "bL33"
    eL33 = "eL33"
    bL34 = "bL34"
    eL34 = "eL34"
    bL35 = "bL35"
    bL36 = "bL36"
    eL36 = "eL36"
    eL37 = "eL37"
    eL38 = "eL38"
    eL39 = "eL39"
    eL40 = "eL40"
    eL41 = "eL41"
    eL42 = "eL42"
    eL43 = "eL43"
    P1P2 = "P1P2"

class MitochondrialRNAClass(str,Enum):
    mtrRNA12S = "mt12SrRNA"  # mitochondrial
    mtrRNA16S = "mt16SrRNA"  # mitochondrial

class CytosolicRNAClass(str,Enum):
    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA"  #  c-bacterial or mitochondrial
    rRNA_23S  = "23SrRNA"  # bacterial
    rRNA_25S  = "25SrRNA"  # plants
    rRNA_5_8S = "5.8SrRNA"  # eukaryotic
    rRNA_18S  = "18SrRNA"  # eukaryotic
    rRNA_28S  = "28SrRNA"  # eukaryotic

class ElongationFactorClass(str,Enum):
    # Eukaryotic
    eEF1A = "eEF1A"
    eEF1B = "eEF1B"
    eFSec = "eFSec"
    eEF2  = "eEF2"
    mtEF4 = "mtEF4"
    eIF5A = "eIF5A"
    eEF3  = "eEF3"
    # Bacterial
    EF_Tu = "EF-Tu"
    EF_Ts = "EF-Ts"
    SelB  = "SelB"
    EF_G  = "EF-G"
    EF4   = "EF4"
    EF_P  = "EF-P"
    Tet_O = "Tet_O"
    Tet_M = "Tet_M"
    RelA  = "RelA"
    BipA  = "BipA"
    # Archaeal
    aEF1A = "aEF1A"
    aEF2  = "aEF2"

class InitiationFactorClass(str,Enum):
    #!Eukaryotic
    eIF1  = "eIF1"
    eIF1A = "eIF1A"

    eIF2_alpha = "eIF2_alpha"
    eIF2_beta  = "eIF2_beta"
    eIF2_gamma = "eIF2_gamma"

    eIF2B_alpha   = "eIF2B_alpha"
    eIF2B_beta    = "eIF2B_beta"
    eIF2B_gamma   = "eIF2B_gamma"
    eIF2B_delta   = "eIF2B_delta"
    eIF2B_epsilon = "eIF2B_epsilon"

    eIF3_subunitA = "eIF3_subunitA"
    eIF3_subunitB = "eIF3_subunitB"
    eIF3_subunitC = "eIF3_subunitC"
    eIF3_subunitD = "eIF3_subunitD"
    eIF3_subunitE = "eIF3_subunitE"
    eIF3_subunitF = "eIF3_subunitF"
    eIF3_subunitG = "eIF3_subunitG"
    eIF3_subunitH = "eIF3_subunitH"
    eIF3_subunitI = "eIF3_subunitI"
    eIF3_subunitJ = "eIF3_subunitJ"
    eIF3_subunitK = "eIF3_subunitK"
    eIF3_subunitL = "eIF3_subunitL"
    eIF3_subunitM = "eIF3_subunitM"

    eIF4F_4A = "eIF4F_4A"
    eIF4F_4G = "eIF4F_4G"
    eIF4F_4E = "eIF4F_4E"

    eIF4B = "eIF4B"
    eIF5B = "eIF5B"
    eIF5  = "eIF5"

    #!Bacterial
    IF1 = "IF1"
    IF2 = "IF2"
    IF3 = "IF3"

    #!Archaeal
    aIF_1A       = "aIF1A"
    aIF_2_alpha  = "aIF2_alpha"
    aIF_2_beta   = "aIF2_beta"
    aIF_2_gamma  = "aIF2_gamma"
    aIF_2B_alpha = "aIF2B_alpha"
    aIF_2B_beta  = "aIF2B_beta"
    aIF_2B_delta = "aIF2B_delta"
    aIF5A        = "aIF5A"
    aIF5B        = "aIF5B"


# LifecycleFactorClass = typing.Union[ElongationFactorClass, InitiationFactorClass]
LifecycleFactorClass = enum_union(ElongationFactorClass, InitiationFactorClass)
ProteinClass         = enum_union(CytosolicProteinClass,MitochondrialProteinClass )
PolypeptideClass     = enum_union(LifecycleFactorClass, ProteinClass)
PolynucleotideClass  = enum_union(CytosolicRNAClass, MitochondrialRNAClass, tRNA)
PolymerClass         = enum_union(PolynucleotideClass, PolypeptideClass)


# ? ----------------------------------------------{ Object Types }------------------------------------------------



class Polymer(BaseModel):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    def to_SeqRecord(self) -> SeqRecord:
        return SeqRecord(
            seq         = Seq(self.entity_poly_seq_one_letter_code_can),
            id          = f"{self.src_organism_ids[0]}",
            description = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id),
            name        = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id)
        )

    assembly_id: int

    asym_ids    : list[str]
    auth_asym_id: str

    parent_rcsb_id: str

    src_organism_names : list[str]
    host_organism_names: list[str]

    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description              : Optional[str] = None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    nomenclature                       : list[PolymerClass]

    @field_serializer('nomenclature')
    def serialize_nomenclature(self, nomenclature_classes: list[PolymerClass], ):
        return [x.value for x in nomenclature_classes]

class Protein(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    @staticmethod
    def from_polymer(p: Polymer, **kwargs):
        if kwargs["pfams"] != None and len(kwargs["pfams"]) > 0:
            pfam_comments = list(
                set([pfam["rcsb_pfam_comment"] for pfam in kwargs["pfams"]])
            )
            pfam_descriptions = list(
                set([pfam["rcsb_pfam_description"] for pfam in kwargs["pfams"]])
            )
            pfam_accessions = list(
                set([pfam["rcsb_pfam_accession"] for pfam in kwargs["pfams"]])
            )

        else:
            pfam_comments = []
            pfam_descriptions = []
            pfam_accessions = []

        return Protein(
            **{
                **p.model_dump(),
                "pfam_accessions": pfam_accessions,
                "pfam_comments": pfam_comments,
                "pfam_descriptions": pfam_descriptions,
                "uniprot_accession": [entry["rcsb_id"] for entry in kwargs["uniprots"]]
                if kwargs["uniprots"] != None and len(kwargs["uniprots"]) > 0
                else [],
            }
        )

    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]
    uniprot_accession: list[str]

    def to_polymer(self) -> Polymer:
        return Polymer(**self.model_dump())

class RNA(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    # pass

class NonpolymericLigand(BaseModel):
    # def metadatum(self) -> NonpolymericLigandMetadatum:
    #     return NonpolymericLigandMetadatum(**self.model_dump())

    class NonpolymerComp(BaseModel):

        class Drugbank(BaseModel):
            class DrugbankInfo(BaseModel):
                cas_number: Optional[str] =None
                description: Optional[str] =None

            class DrugbankContainerIdentifiers(BaseModel):
                drugbank_id: str

            drugbank_container_identifiers: DrugbankContainerIdentifiers
            drugbank_info: DrugbankInfo

        class RcsbChemCompTarget(BaseModel):
            interaction_type: Optional[str] = None
            name: Optional[str] =None
            provenance_source: Optional[str] =None
            reference_database_accession_code: Optional[str] =None
            reference_database_name: Optional[str] =None

        drugbank             : Optional[Drugbank]=None
        rcsb_chem_comp_target: Optional[list[RcsbChemCompTarget]]=None

    chemicalId    : str
    chemicalName  : str
    formula_weight: Optional[float] =None

    pdbx_description   : str
    number_of_instances: int

    nonpolymer_comp: Optional[NonpolymerComp] = None

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
            entity_id: str
            auth_asym_id: str
            auth_seq_id: str

        rcsb_nonpolymer_entity_instance_container_identifiers: NonpolymerEntityInstanceContainerIdentifiers

    class PolymerEntityInstance(BaseModel):
        class PolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str

        rcsb_polymer_entity_instance_container_identifiers: PolymerEntityInstanceContainerIdentifiers

    rcsb_id: str  # ex. 5AFI-1
    nonpolymer_entity_instances: Optional[list[NonpolymerEntityInstance]] =None
    polymer_entity_instances: list[PolymerEntityInstance]

class NomenclatureItem(BaseModel):
    nomenclature: list[str]

class NomenclatureTable(BaseModel):
    __pydantic_root_model__: Dict[str, NomenclatureItem]

class PTCInfo(BaseModel):

    site_9_residues      : list[tuple[str, int]]
    LSU_rRNA_auth_asym_id: str
    midpoint_coordinates : tuple[float, float, float]
    nomenclature_table   : NomenclatureTable


class RibosomeStructure(BaseModel):

    def get_nomenclature_map(self) -> dict[str, list[Polymer]]:
        "Return a map from auth_asym_id to nomenclatures of all polymers in the structure"
        _ = {}

        for prp in self.proteins:
            _[prp.auth_asym_id] = [ x.name for x in  prp.nomenclature]

        for prna in self.rnas:
            _[prna.auth_asym_id] =[ y.name for y in  prna.nomenclature]

        return _

    def __hash__(self):
        return hash(self.rcsb_id)

    rcsb_id   : str
    expMethod : str
    resolution: float

    deposition_date:str

    pdbx_keywords     : Optional[str] =None
    pdbx_keywords_text: Optional[str] = None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year        : Optional[int] = None
    citation_rcsb_authors: Optional[list[str]] = None
    citation_title       : Optional[str] = None
    citation_pdbx_doi    : Optional[str] = None

    # deposition_date    : Optional[int] = None #TODO <-- This is a better source than citation_year

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    assembly_map: Optional[list[AssemblyInstancesMap]] = None
    mitochondrial: bool

    # proteins            : list[Any]
    # rnas                : list[Any]
    # nonpolymeric_ligands: list[Any]
    # polymeric_factors   : list[Any]

    proteins: list[Protein]
    rnas    : list[RNA]

    # polymeric_factors   : list[LifecycleFactor]
    # ? This includes DNA-RNA hybrid strands, DNA and all other polymers
    other_polymers      : list[Polymer]
    nonpolymeric_ligands: list[NonpolymericLigand]

    # TODO: Deprecate this
    # @staticmethod
    # def from_json_profile(d: Any):
    #     return RibosomeStructure(**d)
