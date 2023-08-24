import os
from ribctl.etl.ribosome_assets import RibosomeAssets
from ete3 import NCBITaxa
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl import RIBETL_DATA



def is_descendant_of(taxid: int, struct: str):
   ncbi = NCBITaxa()
   src, hst = RibosomeAssets(struct).get_taxids()
   lineage = ncbi.get_lineage(src[0])
   if lineage is None:
       raise LookupError
   return False if taxid not in lineage else True

def filter_by_parent_tax(taxid:int):
    all_structs = os.listdir(RIBETL_DATA)
    descendants = list(filter(lambda x: is_descendant_of(taxid, x), all_structs))
    return descendants

def __node_lineage(node):
    return NCBITaxa().get_lineage(node.taxid)

def __lift_rank_to_species(taxid: int) -> int:
    """Given a taxid, make sure that it's a SPECIES (as opposed to strain, subspecies, isolate, norank etc.)"""
    ncbi = NCBITaxa()
    if ncbi.get_rank([taxid])[taxid] == 'species':
        return taxid

    else:
        lin = iter(ncbi.get_lineage(taxid))
        node = 1
        while ncbi.get_rank([node])[node] != 'species':
            node = next(lin)
        return node

def tax_classify_struct_proportions(ribosome:RibosomeStructure)->int:

    ids = []
    if ribosome.rnas is not None:
        for rna in ribosome.rnas:
            ids = [*rna.src_organism_ids, *ids]
    
    for protein in ribosome.proteins:
       ids = [*protein.src_organism_ids, *ids]
   
    proportions = {}
    for i in set(ids):
        proportions[i] = ids.count(i)/len(ids)

    return max(proportions, key=proportions.get)