import os
from prody import MSA
from api.ribctl.lib.types.types_ribosome import ProteinClass
from api.ribctl.lib.msalib import msa_profiles_dict_prd
from ete3 import NCBITaxa



def node_lineage(node):
    return NCBITaxa().get_lineage(node.taxid)

def lift_rank_to_species(taxid: int) -> int:
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
