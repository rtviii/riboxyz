from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
import sys
from api.ribxz_api.db_queries import Neo4jQuery
from neo4j_ribosome.adapter import Neo4jBuilder
from ribctl.etl.ribosome_assets import RibosomeAssets
sys.dont_write_bytecode = True
from dotenv import load_dotenv
load_dotenv('.env')



#* Recipe for initializing a new instance from the RIBETL_DATA pool
#* - assumes the RiboosomeStrucutre profiles are rendered



def full_upload():
    adapter = Neo4jBuilder("bolt://localhost:7678", "neo4j", "")
    adapter.initialize_new_instance()
    futures:list[Future] =  []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in RibosomeAssets.list_all_structs()[:10]:
            fut          = executor.submit(partial(adapter.add_structure, rcsb_id, False))
            futures.append(fut)
        wait(futures, return_when=ALL_COMPLETED)

#TODO :  ------- SANITY CHECKS



# with open('./neo4j_adapter/3J7Z.json') as f:
#     rs = RibosomeStructure.model_validate(json.load(f))
#     pprint(rs)


# # print(adapter.see_current_auth())
# # adapter.init_polymer_classes()
# # print(adapter.get_any())
# adapter.add_structure('3j7z')
# adapter.add_structure('4ug0')
# adapter.add_structure('5afi')
# adapter.add_structure('7k00')
# # adapter.sync_with_rcsb(10)
# phn = PhylogenyNode.from_taxid(9605)
# adapter.create_lineage(9606)

# {
#     value: '2',
#     title: 'Bacteria',
#     children: [
#       {
#         value: '66',
#         title: 'parent 1-0',
#         children: [
#           {
#             value: '44',
#             title: 'my leaf',
#           },
#           {
#             value: '22',
#             title: 'your leaf',
#           },
#         ],
#       },
#     ],
#   }

# db = Neo4jQuery()

# global nodes 
# nodes = 0

# def inject_species(node, S:int, F:int, L:int):

#     """This acts on the superkingdom node"""
#     global nodes
#     if node['value'] == F:
#         if len(list(filter(lambda subnode: subnode['value'] == S, node['children'])) ) < 1:
#             node['children'].append({'value': S, 'title': '' })
#             nodes+=1
#     return node

# def inject_families(node, S:int, F:int, K:int):
#     """This acts on the superkingdom node"""
#     global nodes
#     if node['value'] == K:
#         if len(list(filter(lambda subnode: subnode['value'] == F, node['children'])) ) < 1:
#             node['children'].append({'value': F, 'title': [], 'children': []})
#             nodes+=1
#         list(map(lambda node: inject_species(node,S, F, K), node['children']))
#     return node

# for tax in s:
#     p = Taxid.get_lineage(tax, include_only=['superkingdom', 'family', 'species'])
#     K,F,S = p
#     if len( list(filter(lambda obj: obj['value'] == K, normalized_taxa)) ) < 1:
#         normalized_taxa.append({'value': K, 'title': '', "children": []})
#         nodes+=1
    
#     list(map(lambda node: inject_families(node,S, F, K), normalized_taxa))

    

    
# _ = set()
# for i in s:
#     for y in Taxid.get_lineage(i, include_only=['superkingdom', 'family', 'species']):
#         _.add(y)

# print(len(s), nodes, len(_))

    
    # else:
    #     continue

# pprint(normalized_taxa)




# pprint(adapter.init_phylogenies())

