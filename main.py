import json
from pprint import pprint
from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Result, Transaction

from ribctl.lib.types.types_ribosome import RibosomeAssets, RibosomeStructure




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
        # print("\033[91m ------------------------ \033[0m")
        node_lig, node_rib, rel_rel = record["lig"], record["rib"], record['rel']
        # pprint(node_lig.id)
        # pprint(rel_rel.type)
        # pprint(rel_rel.items())
        # pprint(rel_rel.start_node)
        # pprint(rel_rel.end_node)
        # pprint(node_lig.labels)
        # pprint(node_lig.items())
        

#※ ----------------[ 2. Neo4j Data Types ]

# _rib = RibosomeStructure(rcsb_id="4UG0")
# _rib = RibosomeAssets("4UG0").json_profile()
rib  = RibosomeStructure.from_json_profile("4UG0")
r= rib.dict()
# print(json.dumps(RibosomeAssets("4UG0").json_profile()))

# pprint(R)
def create_structure_node(_rib:RibosomeStructure)->Callable[[Transaction|ManagedTransaction], Result]:
    R = _rib.dict()
    def create_parametrized_node(tx:Transaction|ManagedTransaction):
        return  tx.run("""//
    merge ( struct:RibosomeStructure{
              rcsb_id               : $rcsb_id,
              expMethod             : $expMethod,
              resolution            : $resolution,
              citation_rcsb_authors : $citation_rcsb_authors,
              citation_title        : $citation_title,

              pdbx_keywords     : $pdbx_keywords     ,
              pdbx_keywords_text: $pdbx_keywords_text,

              src_organism_ids           : $src_organism_ids,
              src_organism_names         : $src_organism_names,

              host_organism_ids           : $host_organism_ids,
              host_organism_names         : $host_organism_names

              })

              on create set
              struct.rcsb_external_ref_id   = CASE WHEN $rcsb_external_ref_id = null then \"null\" else $rcsb_external_ref_id   END,
              struct.rcsb_external_ref_type = CASE WHEN $rcsb_external_ref_type = null then \"null\" else $rcsb_external_ref_type END,
              struct.rcsb_external_ref_link = CASE WHEN $rcsb_external_ref_link = null then \"null\" else $rcsb_external_ref_link END,
              struct.citation_pdbx_doi      = CASE WHEN $citation_pdbx_doi = null then \"null\" else $citation_pdbx_doi END,
              struct.citation_year          = CASE WHEN $citation_year = null then \"null\" else $citation_year END"
        """,**R)
    return create_parametrized_node
    


                                                            # value.citation_year as cit_year,
                                                            # value.citation_rcsb_authors as cit_authors,
                                                            # value.citation_title as cit_title,
                                                            # value.citation_pdbx_doi as cit_doi,

rib = RibosomeStructure.from_json_profile("4UG0")
r   = rib.dict()


with driver.session() as s:
    s.execute_write(create_structure_node(rib))
