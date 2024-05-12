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

for tax in s:
    p = Taxid.get_lineage(tax, include_only=['superkingdom', 'family', 'species'])
    K,F,S = p
    if len( list(filter(lambda obj: obj['value'] == K, normalized_taxa)) ) < 1:
        print(True)
        normalized_taxa.append({'value': K, 'title': []})
    # else:
    #     continue

print(normalized_taxa)

    



# pprint(adapter.init_phylogenies())

