import os
from prody import MSA
from api.ribctl.lib.types.types_ribosome import ProteinClass
from api.ribctl.msa.msalib import msa_profiles_dict_prd
from ete3 import NCBITaxa
RIBETL_DATA = os.environ.get('RIBETL_DATA')
rcsb_id     = '3J7Z'

ncbi = NCBITaxa()
def get_fasta_taxid(label:str):
    return label.split('|')[-1]

def node_lineage(node):
    return ncbi.get_lineage(node.taxid)

def lift_rank(taxid: int) -> int:
    """Given a taxid, make sure that it's a SPECIES (as opposed to strain, subspecies, isolate, norank etc.)"""
    if ncbi.get_rank([taxid])[taxid] == 'species':
        return taxid
    else:
        lin = iter(ncbi.get_lineage(taxid))
        node = 1
        while ncbi.get_rank([node])[node] != 'species':
            node = next(lin)
        return node

msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()
list_of_taxids = []
bl12 = msa_profiles['bL12']
class_taxids = [get_fasta_taxid(p.getLabel()) for p in bl12]

def msa_pick_taxa(msa: MSA, taxids: list[str]):
    """Given a MSA and a list of taxids, return a new MSA with only the sequences that match the taxids."""
    return MSA([s for s in msa if get_fasta_taxid(s.getLabel()) in taxids])

def phylogenetic_neighborhood(taxids_base:list[str], taxid_target :str, n_neighbors:int=10):
    tree        = ncbi.get_topology([*taxids_base, taxid_target])
    target_node = tree.search_nodes(name=taxid_target)[0]

    phylo_distances   = [(node.name, tree.get_distance(target_node, node)) for node in tree.traverse()]
    phylo_distances_s = sorted(phylo_distances, key=lambda x: x[1])
    nbr_taxids        = list(map(lambda tax_phydist: tax_phydist[0], phylo_distances_s))

    if len(nbr_taxids) < n_neighbors:
        return nbr_taxids
    else:
        return nbr_taxids[:n_neighbors]


# ex.  418699

class_taxids = phylogenetic_neighborhood(class_taxids, '418699', 4)
bl12_picked =  msa_pick_taxa(bl12, class_taxids)
print(bl12)
for m in bl12[:4]:
    print(m)
    print(m.getLabel())

print("--------")
for m in bl12_picked:
    print(m)
    print(m.getLabel())

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