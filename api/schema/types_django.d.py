from ninja import Schema
from ribctl.lib.types.types_ribosome import Protein

class RibosomeResponse(Schema):
    id: int
    name: str
    proteins: list[Protein]
