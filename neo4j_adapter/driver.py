from pprint import pprint
import sys

from api.ribxz_api.db_queries import DBQuery
from ribctl.etl.ribosome_assets import RibosomeAssets
sys.dont_write_bytecode = True
from neo4j_adapter.adapter import Neo4jAdapter
from ribctl.lib.libtax import Taxid, ncbi

from dotenv import load_dotenv

load_dotenv('.env')



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

db = DBQuery()
s = db.get_taxa('source')

import operator
normalized_taxa = []

global nodes 
nodes = 0

def inject_species(node, S:int, F:int, L:int):

    """This acts on the superkingdom node"""
    global nodes
    if node['value'] == F:
        if len(list(filter(lambda subnode: subnode['value'] == S, node['children'])) ) < 1:
            node['children'].append({'value': S, 'title': '' })
            nodes+=1
    return node

def inject_families(node, S:int, F:int, K:int):
    """This acts on the superkingdom node"""
    global nodes
    if node['value'] == K:
        if len(list(filter(lambda subnode: subnode['value'] == F, node['children'])) ) < 1:
            node['children'].append({'value': F, 'title': [], 'children': []})
            nodes+=1
        list(map(lambda node: inject_species(node,S, F, K), node['children']))
    return node

for tax in s:
    p = Taxid.get_lineage(tax, include_only=['superkingdom', 'family', 'species'])
    K,F,S = p
    if len( list(filter(lambda obj: obj['value'] == K, normalized_taxa)) ) < 1:
        normalized_taxa.append({'value': K, 'title': '', "children": []})
        nodes+=1
    
    list(map(lambda node: inject_families(node,S, F, K), normalized_taxa))

    

    
_ = set()
for i in s:
    for y in Taxid.get_lineage(i, include_only=['superkingdom', 'family', 'species']):
        _.add(y)

print(len(s), nodes, len(_))

    
    # else:
    #     continue

# pprint(normalized_taxa)

    



# pprint(adapter.init_phylogenies())

