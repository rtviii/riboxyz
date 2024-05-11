import json
from pprint import pprint
from typing import Any, Optional
from enum import Enum
import typing
from typing_extensions import Literal
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pydantic import BaseModel, field_serializer
from ribctl import NCBI_TAXA_SQLITE
from ribctl.lib.enumunion import enum_union
from ete3 import NCBITaxa

# TODO:
# |********************************************************************************************************|
# | https://docs.google.com/spreadsheets/d/1mapshbn1ArofPN-Omu8GG5QdcwlJ0ym0BlN252kkUBU/edit#gid=815712128 |
# |********************************************************************************************************|
# ? ---------------- [ Phylogeny] ----------------

TAXID_BACTERIA  = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA   = 2157
PhylogenyRank = Literal["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
ncbi          = NCBITaxa(dbfile=NCBI_TAXA_SQLITE)
class Taxid:
    @staticmethod
    def is_descendant_of(parent_taxid: int, target_taxid: int) -> bool:
        lineage = ncbi.get_lineage(target_taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        return False if parent_taxid not in lineage else True

    @staticmethod
    def get_name(taxid):
        return list(ncbi.get_taxid_translator([taxid]).values())[0]

    @staticmethod
    def get_lineage(taxid)->list[int]:
        """Return ncbi lineage, except filter out the ranks that are not among the @PhylogenyRank."""
        print("Getting lineage for", ncbi.get_rank(ncbi.get_lineage(taxid)))
        pprint(ncbi.get_rank(list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )))
        return list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )


    @staticmethod
    def rank(taxid: int) -> PhylogenyRank:
        """Given a @taxid, return the rank of the taxid"""
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage)[taxid]

    @staticmethod
    def coerce_to_rank(taxid: int, target_rank: PhylogenyRank) -> int | None:
        """Given a @taxid and a @rank, return the taxid of the first ancestor of @taxid that is at @rank"""
        lineage = ncbi.get_lineage(taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        for item in lineage:
            rank = ncbi.get_rank([item])[item]
            if rank == target_rank:
                return item

        raise IndexError("Taxid {} does not have a {} level".format(taxid, target_rank))

    @staticmethod
    def superkingdom(
        taxid: int,
    ) -> typing.Literal["bacteria", "eukaryota", "archaea", "virus"]:
        match (
            Taxid.is_descendant_of(TAXID_EUKARYOTA, taxid),
            Taxid.is_descendant_of(TAXID_BACTERIA, taxid),
            Taxid.is_descendant_of(TAXID_ARCHAEA, taxid),
        ):
            case (False, False, True):
                return "archaea"
            case (False, True, False):
                return "bacteria"
            case (True, False, False):
                return "eukaryota"
            case (False, False, False):
                print("Probably a virus")
                return "virus"
            case _:
                raise ValueError(
                    "Taxid {} is not a descendant of any of the three domains".format(
                        taxid
                    )
                )

    @staticmethod
    def coerce_all_to_rank(taxids: list[int], level: PhylogenyRank) -> list[int]:
        """Given a list of taxids, return a list of the same taxids but coerced to the given lineage level(rank)."""
        new = []
        for taxid in taxids:
            try:
                new.append(Taxid.coerce_to_rank(taxid, level))
            except Exception as e:
                print(e)
                raise Exception(e)
        return new

    @staticmethod
    def get_descendants_of(parent: int, targets: list[int]):
        """Given a @parent taxid and a list of @taxids, return the subset of @taxids that are descendants of @parent"""
        descendants = set()
        for tax_id in targets:
            lineage = ncbi.get_lineage(tax_id)
            if lineage == None:
                raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
            if parent in lineage:
                descendants.add(tax_id)
        return descendants

class PhylogenyNode(BaseModel):

    @staticmethod
    def from_taxid(taxid:int):
        return PhylogenyNode(ncbi_tax_id=taxid, scientific_name=Taxid.get_name(taxid), rank=Taxid.rank(taxid))

    def __hash__(self) -> int:
        return self.ncbi_tax_id

    def get_lineage(self) -> list[int]:
        return Taxid.get_lineage(self.ncbi_tax_id)

    def get_lineage_nodes(self) -> list[int]:
        _ = []
        for taxid in self.get_lineage():
            if Taxid.rank(taxid) not in typing.get_args(PhylogenyRank):
               continue
            _.append(self.from_taxid(taxid))
        return _

    ncbi_tax_id    : int
    scientific_name: str
    rank           : PhylogenyRank


# ? ----------------------------------------------{ Subcomponent Types }------------------------------------------------
class tRNA(Enum):
    tRNA = "tRNA"

class MitochondrialProteinClass(Enum):
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

class CytosolicProteinClass(Enum):
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

class MitochondrialRNAClass(Enum):
    mtrRNA12S = "mt12SrRNA"  # mitochondrial
    mtrRNA16S = "mt16SrRNA"  # mitochondrial

class CytosolicRNAClass(Enum):
    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA"  #  c-bacterial or mitochondrial
    rRNA_23S  = "23SrRNA"  # bacterial
    rRNA_25S  = "25SrRNA"  # plants
    rRNA_5_8S = "5.8SrRNA"  # eukaryotic
    rRNA_18S  = "18SrRNA"  # eukaryotic
    rRNA_28S  = "28SrRNA"  # eukaryotic

class ElongationFactorClass(Enum):
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

class InitiationFactorClass(Enum):
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

# class PolymerMetadatum(BaseModel):
#     assembly_id           : int
#     asym_ids              : list[str]
#     auth_asym_id          : str
#     entity_poly_seq_length: int
#     nomenclature          : list[PolymerClass]

#     @field_serializer('nomenclature')
#     def serialize_dt(self, ncl: list[PolymerClass], _info):
#         return [x.value for x in ncl]

class Polymer(BaseModel):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    def enum_union_fix(_):
        return _
    def to_dict(self):
        """A hack for enum.union to work with pydantic BaseModel. Otherwise EnumUnion instances are represented as <MitochondrialProteinClass.mL64: 'mL64'> etc.(Correct is "mL64")"""
        return json.loads(self.model_dump_json())

    def to_SeqRecord(self) -> SeqRecord:
        return SeqRecord(
            seq         = Seq(self.entity_poly_seq_one_letter_code_can),
            id          = f"{self.src_organism_ids[0]}",
            description = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id),
            name        = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id)
        )

    # def metadatum(self) -> PolymerMetadatum:

    #     return PolymerMetadatum(
    #         assembly_id            = self.assembly_id,
    #         asym_ids               = self.asym_ids,
    #         auth_asym_id           = self.auth_asym_id,
    #         entity_poly_seq_length = self.entity_poly_seq_length,
    #         nomenclature           = [x.value for x in self.nomenclature])

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
    def serialize_dt(self, ncl: list[PolymerClass], _info):
        return [x.value for x in ncl]

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
        class ChemComp(BaseModel):
            class ChemCompContainerIdentifiers(BaseModel):
                id: str
                name: str
                three_letter_code: str

            chem_comp_container_identifiers: ChemCompContainerIdentifiers

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

        chemp_comp           : Optional[ChemComp]=None
        drugbank             : Optional[Drugbank]=None
        rcsb_chem_comp_target: Optional[list[RcsbChemCompTarget]]=None

    chemicalId    : str
    chemicalName  : str
    formula_weight: Optional[float] =None

    pdbx_description: str
    number_of_instances: int

    nonpolymer_comp: Optional[NonpolymerComp] =None


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

# class RibosomeStructureMetadatum(BaseModel):

#     rcsb_id   : str
#     expMethod : str
#     resolution: float

#     pdbx_keywords     : Optional[str] =None
#     pdbx_keywords_text: Optional[str] = None

#     rcsb_external_ref_id: list[str]
#     rcsb_external_ref_type: list[str]
#     rcsb_external_ref_link: list[str]

#     citation_year         : Optional[int]      = None
#     citation_rcsb_authors : Optional[list[str]] = None
#     citation_title        : Optional[str]      = None
#     citation_pdbx_doi     : Optional[str]      = None

#     src_organism_ids  : list[int]
#     src_organism_names: list[str]

#     host_organism_ids  : list[int]
#     host_organism_names: list[str]

#     # assembly_map: list[AssemblyInstancesMap]
#     mitochondrial: bool

#     # proteins_metadata: list[PolymerMetadatum]
#     # rnas_metadata    : list[PolymerMetadatum]
#     # # ? This includes DNA-RNA hybrid strands, DNA and all other polymers
#     # other_polymers_metadata: list[PolymerMetadatum]
#     # nonpolymeric_ligands   : list[NonpolymericLigandMetadatum]

class RibosomeStructure(BaseModel):

    def get_nomenclature_map(self) -> dict[str, list[Polymer]]:
        "Return a map from auth_asym_id to nomenclatures of all polymers in the structure"
        _ = {}

        for prp in self.proteins:
            _[prp.auth_asym_id] = [ x.name for x in  prp.nomenclature]

        for prna in self.rnas:
            _[prna.auth_asym_id] =[ y.name for y in  prna.nomenclature]

        return _


    rcsb_id   : str
    expMethod : str
    resolution: float

    pdbx_keywords     : Optional[str] =None
    pdbx_keywords_text: Optional[str] = None

    rcsb_external_ref_id: list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year        : Optional[int] = None
    citation_rcsb_authors: Optional[list[str]] = None
    citation_title       : Optional[str] = None
    citation_pdbx_doi    : Optional[str] = None

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids: list[int]
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
