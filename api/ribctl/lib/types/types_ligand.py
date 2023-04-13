from pydantic import BaseModel

from api.ribctl.lib.types.types_ribosome import Polymer



class LigandLikePolymer(Polymer):


class Ligand(BaseModel)  : 

      chemicalId         : str
      chemicalName       : str
      formula_weight     : None | float
      pdbx_description   : str
      number_of_instances: int
