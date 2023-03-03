from pprint import pprint
from ribctl.lib.types.types_ribosome import RibosomeAssets, RibosomeStructure


# To generate types:
# pydantic2ts --module ribctl/lib/types/types_ribosome.py --output ./types_ribosome.ts
r = RibosomeStructure.from_json_profile("4UG0")
pprint(r.rnas[0].dict())
      