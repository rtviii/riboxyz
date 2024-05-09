from neo4j_adapter.adapter import Neo4jAdapter
from dotenv import load_dotenv

load_dotenv('.env')


adapter = Neo4jAdapter('bolt://localhost:7687', 'neo4j')
# print(adapter.see_current_auth())
# adapter.init_polymer_classes()

# print(adapter.get_any())
adapter.add_structure('3j7z')
    
# adapter.sync_with_rcsb(10)


