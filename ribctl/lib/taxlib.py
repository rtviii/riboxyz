from typing import Tuple
import typing
from ete3 import NCBITaxa
from ribctl import RIBETL_DATA, TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_EUKARYOTA


"""
These methods are intended to work across structures and their components. 
Primary operations should be in terms of integer tax. ids not objects.

Separately implement the source/host thing for structs.
"""






def taxid_superkingdom(
    taxid: int,
) -> typing.Literal["bacteria", "eukaryota", "archaea", "virus"]:
    match (
        taxid_is_descendant_of(TAXID_EUKARYOTA, taxid)[0],
        taxid_is_descendant_of(TAXID_BACTERIA, taxid)[0],
        taxid_is_descendant_of(TAXID_ARCHAEA, taxid)[0],
    ):
        case (False, False, True):
            return "archaea"
        case (False, True, False):
            return "bacteria"
        case (True, False, False):
            return "eukaryota"
        case (False, False, False):
            print("Probably a virus")
            return "virus"
        case _:
            raise ValueError(
                "Taxid {} is not a descendant of any of the three domains".format(taxid)
            )


def get_descendants_of(parent: int, targets: list[int]):
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

def taxid_is_descendant_of(
    parent_taxid: int, target_taxid: int
) -> (bool, list[int] | None):
    ncbi = NCBITaxa()
    lineage = ncbi.get_lineage(target_taxid)
    if lineage is None:
        raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
    return (False, lineage) if parent_taxid not in lineage else (True, lineage)


def descendants_of_taxid(
    struct_taxids: list[Tuple[str, int]], parent_taxid: int
) -> list[Tuple[str, int]]:
    descendants = []
    for rcsb_id, taxid in struct_taxids:
        descends, lin = taxid_is_descendant_of(parent_taxid, taxid)
        if descends:
            descendants.append((rcsb_id, taxid))

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