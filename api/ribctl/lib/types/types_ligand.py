from typing import Literal
from pydantic import BaseModel

from api.ribctl.lib.types.types_ribosome import Polymer

LigandlikePolymerClass = Literal[
      "Elongation Factor",
      "Initiation Factor",
      "Translation Factor"
]


class LigandLikePolymer(Polymer):
      nomenclature:
      




class Ligand(BaseModel)  : 

      chemicalId         : str
      chemicalName       : str
      formula_weight     : None | float
      pdbx_description   : str
      number_of_instances: int
