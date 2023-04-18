from pydantic import BaseModel
from .types_poly_nonpoly_ligand import LSU_Proteins, NonpolymericLigandClass, PolymericFactorClass, RNAClass, SSU_Proteins


class Organism(BaseModel): 
      domain             : tuple[int, str]
      species            : tuple[int, str]
      strain             : tuple[int, str]




