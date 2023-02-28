from pprint import pprint
from ribctl.lib.types.types_ribosome import RibosomeAssets, RibosomeStructure



r = RibosomeStructure.from_json_profile("4UG0")
pprint(r.rnas[0].dict())
      