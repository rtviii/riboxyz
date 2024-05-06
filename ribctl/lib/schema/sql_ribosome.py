from typing import Optional
from sqlalchemy import create_engine
from sqlmodel import Field, SQLModel, Session
import sys
sys.path.append('/home/rtviii/dev/riboxyz/')
from ribctl.lib.schema.types_ribosome import PolymerClass, PolynucleotideClass


class RibosomeStructureMetadatum(SQLModel, table=True):

    rcsb_id   : str  = Field(primary_key=True)
    expMethod : str
    resolution: float

    pdbx_keywords     : Optional[str] =None
    pdbx_keywords_text: Optional[str] = None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year         : Optional[int]      = None
    citation_rcsb_authors: Optional[list[str]] = None
    citation_title        : Optional[str]      = None
    citation_pdbx_doi     : Optional[str]      = None

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    mitochondrial: bool

class Polymer(SQLModel, table=True):

    parent_rcsb_id: str = Field(foreign_key="RibosomeStructureMetadatum.rcsb_id")
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
    nomenclature                       : list[str]



