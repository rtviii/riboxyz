import os
from ribctl.etl.ribosome_assets import RibosomeAssets
from ete3 import NCBITaxa
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl import RIBETL_DATA


"""
These methods are intended to work across structures and their components. 
Primary operations should be in terms of integer tax. ids not objects.

Separately implement the source/host thing for structs.
"""




def get_descendants_of( parent:int, targets:list[int]):
    ncbi = NCBITaxa()
    
    # Get the taxonomic lineage of the parent tax id
    parent_lineage = ncbi.get_lineage(parent)
    
    # Get the descendants of the parent tax id
    descendants = set()

    for tax_id in targets:
        lineage = ncbi.get_lineage(tax_id)
        if parent in lineage:
            descendants.add(tax_id)
    
    print(descendants)
    return descendants

# ? Struct-specific functions

def is_descendant_of(taxid: int, struct: str) -> bool:
    ncbi = NCBITaxa()
    src, hst = RibosomeAssets(struct).get_taxids()
    lineage = ncbi.get_lineage(src[0])
    if lineage is None:
        raise LookupError
    return False if taxid not in lineage else True


def filter_by_parent_tax(taxid: int):
    all_structs = os.listdir(RIBETL_DATA)
    descendants = list(filter(lambda x: is_descendant_of(taxid, x), all_structs))
    return descendants


def __node_lineage(node):
    return NCBITaxa().get_lineage(node.taxid)


def __lift_rank_to_species(taxid: int) -> int:
    """Given a taxid, make sure that it's a SPECIES (as opposed to strain, subspecies, isolate, norank etc.)"""
    ncbi = NCBITaxa()
    if ncbi.get_rank([taxid])[taxid] == "species":
        return taxid

    else:
        lin = iter(ncbi.get_lineage(taxid))
        node = 1
        while ncbi.get_rank([node])[node] != "species":
            node = next(lin)
        return node


def classify_struct_by_proportions(ribosome: RibosomeStructure) -> int:
    ids = []
    if ribosome.rnas is not None:
        for rna in ribosome.rnas:
            ids = [*rna.src_organism_ids, *ids]

    for protein in ribosome.proteins:
        ids = [*protein.src_organism_ids, *ids]

    proportions = {}
    for i in set(ids):
        proportions[i] = ids.count(i) / len(ids)

    return max(proportions, key=proportions.get)
