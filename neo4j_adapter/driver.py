import sys
sys.dont_write_bytecode = True
import json
from pprint import pprint
from neo4j_adapter.adapter import Neo4jAdapter
from ribctl.lib.schema.types_ribosome import RibosomeStructure

from dotenv import load_dotenv

load_dotenv('.env')



# with open('./neo4j_adapter/3J7Z.json') as f:
#     rs = RibosomeStructure.model_validate(json.load(f))
#     pprint(rs)


adapter = Neo4jAdapter('bolt://localhost:7687', 'neo4j')

# # print(adapter.see_current_auth())
# # adapter.init_polymer_classes()
# # print(adapter.get_any())
adapter.add_structure('3j7z')
adapter.add_structure('4ug0')
adapter.add_structure('5afi')
adapter.add_structure('7k00')
# # adapter.sync_with_rcsb(10)


