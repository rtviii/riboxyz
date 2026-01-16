from typing import Dict, Optional
from enum import Enum
import typing
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.PDB.Residue import Residue
from pydantic import BaseModel
from typing import Optional
from ribctl.lib.types.polymer import PolymerClass
from ribctl.lib.schema.primitives import AMINO_ACIDS, NUCLEOTIDES


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

class Protein(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    @staticmethod
    def from_polymer(p: Polymer, **kwargs):
        if kwargs["pfams"] != None and len(kwargs["pfams"]) > 0:

            pfam_comments     = list( set([pfam["rcsb_pfam_comment"] for pfam in kwargs["pfams"]]) )
            pfam_descriptions = list( set([pfam["rcsb_pfam_description"] for pfam in kwargs["pfams"]]) )
            pfam_accessions   = list( set([pfam["rcsb_pfam_accession"] for pfam in kwargs["pfams"]]) )

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
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
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

    SMILES       : Optional[str] =None
    SMILES_stereo: Optional[str] =None
    InChI        : Optional[str] =None
    InChIKey     : Optional[str] =None

class ResidueSummary(BaseModel): 
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }

    label_seq_id : typing.Optional[int] = None
    label_comp_id: typing.Optional[str] = None
    auth_asym_id : str
    auth_seq_id  : int
    rcsb_id      : str
    full_id      : typing.Optional[tuple[str, int, str, tuple[str, int, str]]]


    @staticmethod
    def is_canonical(resname:str):
        return resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]

    @staticmethod
    def three_letter_code_to_one(resname: str):
        if resname in AMINO_ACIDS:
            return AMINO_ACIDS[resname]["one_letter_code"]
        elif resname in NUCLEOTIDES:
                return resname
        else:
            return '-'

    @staticmethod
    def one_letter_code_to_three(resname: str):
        if resname in [*map(lambda x: x[1]['one_letter_code'], AMINO_ACIDS.items())]:
            for tlk, d in AMINO_ACIDS.items():
                if d["one_letter_code"] == resname:
                    return tlk
        elif resname in NUCLEOTIDES:
                return resname
        else:
            return '-'

    def __hash__(self):
        return hash( self.get_resname() if self.get_resname() is not None else "" + str(self.get_seqid()) + self.get_parent_auth_asym_id() )

    def get_resname(self):
        return self.label_comp_id

    def get_seqid(self):
        (structure_id, model_id, chain_id, _) = self.full_id
        (hetero, seqid, insertion_code)       = _
        return seqid

    def get_parent_auth_asym_id(self):
        (structure_id, assembly_id, chain_id, _) = self.full_id
        return chain_id

    @staticmethod
    def from_biopython_residue(r: Residue):

        (structure_id, model_id, chain_id, _) = r.get_full_id()
        (hetero, seqid, insertion_code) = _

        return ResidueSummary(
            auth_seq_id   = seqid,
            label_seq_id  = None,
            label_comp_id = r.get_resname(),
            auth_asym_id  = chain_id,
            full_id       = r.get_full_id(),
            rcsb_id       = r.get_full_id()[0]
        )

class AssemblyInstancesMap(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
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

    rcsb_id                    : str  # ex. 5AFI-1
    nonpolymer_entity_instances: Optional[list[NonpolymerEntityInstance]] =None
    polymer_entity_instances   : list[PolymerEntityInstance]

class NomenclatureItem(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    nomenclature: list[str]

class NomenclatureTable(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    __pydantic_root_model__: Dict[str, NomenclatureItem]

class PTCInfo(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    location: list[float]
    residues: list[ResidueSummary]

class ConstrictionSite(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    location: list[float]

class RibosomeStructureMetadata(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }

    def get_polymers_by_assembly(self)->dict[str,list[str]]:
        if self.assembly_map:
            _ = {}
            for assembly in self.assembly_map:
                _[assembly.rcsb_id] = [poly.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id for poly in assembly.polymer_entity_instances]
            return _
        else:
            raise Exception("No assembly map found")

    def get_nomenclature_map(self) -> dict[str, list[PolymerClass]]:
        "Return a map from auth_asym_id to nomenclatures of all polymers in the structure"
        _ = {}

        for prp in self.proteins:
            _[prp.auth_asym_id] = [ x.name for x in  prp.nomenclature]

        for prna in self.rnas:
            _[prna.auth_asym_id] =[ y.name for y in  prna.nomenclature]

        return _




    def __hash__(self):
        return hash(self.rcsb_id)

    @staticmethod
    def model_validate_partial(data: dict, polymer_fields: bool = False):
        if polymer_fields:
            return RibosomeStructure.model_validate(data)
        else:
            # Create a copy of the data without the detailed fields
            filtered_data = {k: v for k, v in data.items() if k not in ['proteins', 'rnas', 'other_polymers', 'nonpolymeric_ligands']}
            return RibosomeStructure.model_validate(filtered_data)

    rcsb_id   : str
    expMethod : str
    resolution: float

    deposition_date:Optional[ str ] = None

    pdbx_keywords      : Optional[str] = None
    pdbx_keywords_text: Optional[str]  = None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year         : None| Optional[int|str] = None
    citation_rcsb_authors: None|Optional[list[str]] = None
    citation_title        :None| Optional[str]      = None
    citation_pdbx_doi     :None| Optional[str]      = None

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    assembly_map    : Optional[list[AssemblyInstancesMap]]  = None
    mitochondrial   : bool
    subunit_presence: Optional[list[typing.Literal['ssu','lsu']]] = None

class RibosomeStructure(RibosomeStructureMetadata):

    rnas                : list[RNA]
    other_polymers      : list[Polymer]
    nonpolymeric_ligands: list[NonpolymericLigand]
    proteins            : list[Protein]

