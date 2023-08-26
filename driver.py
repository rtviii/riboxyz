"""this is just a place to experiment with things and not pollute the cli."""
from ribctl.lib.util_taxonomy import tax_list_filterby_taxid
from ribctl.ribosome_assets import RibosomeAssets


def test():
    all = RibosomeAssets.list_all_structs()
    s   = []
    for struct in all[:100]:
        s=  [*s , RibosomeAssets(struct).get_taxids()]
    print(s)
    # tax_list_filterby_taxid(all, 2759)
    
