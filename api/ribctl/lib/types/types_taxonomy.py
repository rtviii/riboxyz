from typing import Any
import typing
from pydantic import BaseModel
from .types_poly_nonpoly_ligand import LSU_Proteins, NonpolymericLigandClass, PolymericFactorClass, RNAClass, SSU_Proteins


class Organism(BaseModel):
    asym_ids: list[str]
    auth_asym_id: str

    parent_rcsb_id   : str

    src_organism_names : list[str]
    host_organism_names: list[str]
    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description: str | None