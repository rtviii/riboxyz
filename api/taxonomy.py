import os
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import ProteinClass
from api.ribctl.lib.msalib import msa_profiles_dict_prd
from ete3 import NCBITaxa
from api.scratch_tunnel_workflow import RIBETL_DATA

ncbi = NCBITaxa()

def is_descendant_of(taxid: int, struct: str):
   src, hst = RibosomeAssets(struct).get_taxids()
   lineage = ncbi.get_lineage(src[0])
   if lineage is None:
       raise LookupError
   return False if taxid not in lineage else True

def filter_by_parent_tax(taxid:int):
    all_structs = os.listdir(RIBETL_DATA)
    descendants  = list(filter(lambda x: is_descendant_of(taxid, x), all_structs))
    return descendants

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
