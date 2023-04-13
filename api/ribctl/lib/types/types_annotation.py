from typing import Any
import typing
from pydantic import BaseModel
from .types_polymer_nonpoly_ligand import LSU_Proteins, RNAClass, SSU_Proteins

class InterProFamily(BaseModel):
    family_id: str
    type: str
    description: str

class GOClass(BaseModel):
    class_id: str
    annotation: str

class PFAMFamily(BaseModel):
    family_id: str
    annotation: str
    family_type: str

class LastUpdate(BaseModel):
    date: str
    added_structure: str