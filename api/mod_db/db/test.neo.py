from api.mod_db.db.ribosomexyz import ribosomexyzDB
from os import environ

from ribctl.etl.ribosome_assets import RibosomeAssets


print(environ)


rbxz = ribosomexyzDB('bolt://localhost:7687', 'neo4j', 'password')
ra = RibosomeAssets('3j7z')
rbxz.add_structure(ra)
print(rbxz.get_any())


