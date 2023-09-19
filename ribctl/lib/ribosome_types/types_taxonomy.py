from pydantic import BaseModel
from .types_poly_nonpoly_ligand import LSUProteinClass, NonpolymericLigandClass, PolymericFactorClass, RNAClass, SSUProteinClass


class Organism(BaseModel): 
      domain             : tuple[int, str]
      species            : tuple[int, str]
      strain             : tuple[int, str]




