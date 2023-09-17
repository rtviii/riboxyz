"""this is just a place to experiment with things and not pollute the cli."""
from ribctl.lib.util_taxonomy import  get_descendants_of
from ribctl.ribosome_assets import RibosomeAssets
from ete3 import NCBITaxa


def test():
    all = RibosomeAssets.list_all_structs()
    s   = []
    for struct in all[:100]:
        s=  [*s , RibosomeAssets(struct).get_taxids()]
    sourceids = []
    for tup in s:
        sourceids +=tup[0]

    print(sourceids)
    
    #! -------------
    parent = 2759
    targets = sourceids
    ncbi = NCBITaxa()
    
    # Get the taxonomic lineage of the parent tax id
    parent_lineage = ncbi.get_lineage(parent)
    
    # Get the descendants of the parent tax id
    descendants = set()
    for tax_id in targets:
        lineage = ncbi.get_lineage(tax_id)
        if parent in lineage:
            descendants.add(tax_id)
    
    print(get_descendants_of(2759, sourceids))
    return descendants
    # tax_list_filterby_taxid(all, 2759)
    