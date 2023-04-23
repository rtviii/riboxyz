import os
from ete3 import NCBITaxa, NodeStyle, TreeStyle
from api.ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets, obtain_assets, obtain_assets_threadpool

# Create an instance of the NCBITaxa class
# ncbi = NCBITaxa()

RIBETL_DATA = os.environ.get('RIBETL_DATA')

# for struct in os.listdir(RIBETL_DATA):
#     print(RibosomeAssets(struct).profile())

rcsb_id = '3J7Z'


#TODO: Ignore intermediate taxa (induct the species, strains, subspecies etc. as nodes)
# - orgnaism
# Get the taxonomic IDs of all the taxa in the database
# all_taxa =ncbi.get_descendant_taxa(1)
# all_taxa = ncbi.get_descendant_taxa(ncbi.get_descendant_taxa(1))
# print()

# Build a tree object representing the taxonomic hierarchy
# tree = ncbi.get_topology(all_taxa)

# # Customize the appearance of the tree nodes
# ns = NodeStyle()
# ns["size"] = 0

# # Create a TreeStyle object to specify the overall layout and appearance of the tree
# ts = TreeStyle()
# ts.layout_fn = "layout"
# ts.mode = "c"

# # Render the tree to a file
# tree.render("ncbi_taxonomy.png", tree_style=ts)