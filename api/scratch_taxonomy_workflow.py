import os
from pprint import pprint
from prody import MSA
from api.ribctl.lib.types.types_ribosome import ProteinClass
from api.ribctl.msa.msalib import msa_profiles_dict_prd
from ete3 import NCBITaxa
RIBETL_DATA = os.environ.get('RIBETL_DATA')
rcsb_id = '3J7Z'
ncbi = NCBITaxa()

msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()


def get_fasta_taxid(label: str):
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

def msa_pick_taxa(msa: MSA, taxids: list[str])->MSA:
    """Given a MSA and a list of taxids, return a new MSA with only the sequences that match the taxids."""
    seqlabel_tups = [(s, s.getLabel())
                     for s in msa if get_fasta_taxid(s.getLabel()) in taxids]
    seqs = [tup[0] for tup in seqlabel_tups]
    labels = [tup[1] for tup in seqlabel_tups]

    return MSA(seqs, labels=labels)

def phylogenetic_neighborhood(taxids_base: list[str], taxid_target: str, n_neighbors: int = 10)->list[str]:
    tree            = ncbi.get_topology([*taxids_base, taxid_target])
    target_node     = tree.search_nodes(name=taxid_target)[0]
    phylo_all_nodes = [
        (node.name, tree.get_distance(target_node, node))
        for node in tree.traverse()]

    phylo_extant_nodes = filter(
        lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)
    phylo_sorted_nodes = sorted(phylo_extant_nodes, key=lambda x: x[1])
    nbr_taxids = list(
        map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

    # the first element is the target node
    if len(nbr_taxids) < n_neighbors:
        return nbr_taxids[1:] 
    else:
        return nbr_taxids[1:n_neighbors+1]

def msa_yield_taxa(msa: MSA)->list[str]:
    return [get_fasta_taxid(p.getLabel()) for p in msa]

bl12           = msa_profiles['bL12']
class_taxids   = msa_yield_taxa(bl12)
neigbors       = phylogenetic_neighborhood(class_taxids, '418699', 4)
bl12_picked    = msa_pick_taxa(bl12, neigbors)




# TODO: Ignore intermediate taxa (induct the species, strains, subspecies etc. as nodes)
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
