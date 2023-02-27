from pprint import pprint
from neo4j import GraphDatabase, Driver, Result, Transaction


from ribctl.types.types_ribosome import Ligand, RibosomeAssets, RibosomeStructure


def init_driver(uri, username, password):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver("neo4j://localhost:7687",auth=(username, password))
    return api




driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")


#※ -------------------------------= [ 1. Create a Node ]

# # Unit of work
# def get_structures(tx, rcsb_id): # (1)
#     result:Result= tx.run("MATCH (ribosome:RibosomeStructure {rcsb_id: $rcsb_id})-[]-(protein:Protein) RETURN ribosome,protein", rcsb_id=rcsb_id)
#     pprint(result.values('ribosome','protein'))
#     return result

# # Open a Session
# with driver.session() as session:
#     # Run the unit of work within a Read Transaction
#     session.execute_read(get_structures, rcsb_id="5AFI") 
#     session.close()

#※ ----------------[ 2. Neo4j Data Types ]
#### Exploring Records

# When accessing a record, either within a loop, list comprehension or within a single record, you can use the [] bracket syntax.

# The following example extracts the p value fro


with driver.session() as s:
    res = s.run("MATCH (lig:Ligand)-[rel]-(rib:RibosomeStructure) RETURN lig,rib,rel limit 5")
    # for record in res.values():
    #     print("\033[91m ------------------------ \033[0m")
    #     pprint(p)

    for record in res:
        print("\033[91m ------------------------ \033[0m")
        node_lig, node_rib, rel_rel = record["lig"], record["rib"], record['rel']
        # pprint(node_lig.id)
        pprint(rel_rel.type)
        pprint(rel_rel.items())
        pprint(rel_rel.start_node)
        pprint(rel_rel.end_node)
        # pprint(node_lig.labels)
        # pprint(node_lig.items())
        

#※ ----------------[ 2. Neo4j Data Types ]

# _rib = RibosomeStructure(rcsb_id="4UG0")
_rib = RibosomeAssets("4UG0").json_profile()

RibosomeStructure(**_rib)
pprint(_rib)
# pprint(R)
# def create_structure_node(tx:Transaction,_rib:RibosomeStructure):
#     R = _rib.dict()
#     pprint(R)
    # return  tx.run("""""",**)
    




# with driver.session() as s:
#     s.execute_write()
