from pprint import pprint
from typing import List, TypeVar
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
import numpy as np
from ribctl.etl.etl_assets_ops import RibosomeOps
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection

from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import BiopythonChain
from ribctl.lib.schema.types_binding_site import ResidueSummary


rcsb_ids = [
      "2FTC",
      "3J6B",
      "3J7Y",
      "3J9M",
      "3JD5",
      "4CE4",
      "4V19",
      "5AJ3",
      "5AJ4",
      "5MRC",
      "5MRE",
      "5MRF",
      "5OOL",
      "5OOM",
      "6GAW",
      "6GAZ",
      "6GB2",
      "6I9R",
      "6NU2",
      "6NU3",
      "6RW4",
      "6RW5",
      "6VLZ",
      "6VMI",
      "6XYW",
      "6ZM5",
      "6ZM6",
      "7A5F",
      "7A5G",
      "7A5H",
      "7A5I",
      "7A5J",
      "7A5K",
      "7L08",
      "7L20",
      "7O9K",
      "7O9M",
      "7ODR",
      "7ODS",
      "7ODT",
      "7OF0",
      "7OF2",
      "7OF3",
      "7OF4",
      "7OF5",
      "7OF6",
      "7OF7",
      "7OI6",
      "7OI7",
      "7OI8",
      "7OI9",
      "7OIA",
      "7OIB",
      "7OIC",
      "7OID",
      "7OIE",
      "7P2E",
      "7PD3",
      "7PKT",
      "7PNT",
      "7PNU",
      "7PNV",
      "7PNW",
      "7PNX",
      "7PNY",
      "7PNZ",
      "7PO0",
      "7PO1",
      "7PO2",
      "7PO3",
      "7PO4",
      "7QH6",
      "7QH7",
      "7QI4",
      "7QI5",
      "7QI6",
      "8A22",
      "8ANY",
      "8APN",
      "8APO",
      "8CSP",
      "8CSQ",
      "8CSR",
      "8CSS",
      "8CST",
      "8CSU",
      "8OIN",
      "8OIP",
      "8OIQ",
      "8OIR",
      "8OIS",
      "8OIT",
      "8PK0",
      "8QSJ"
    ]


REFERENCE_MITO_STRUCTURE_TRNA             = ( '7A5F' , '24')
# TODO:
# 1.combine the cterm/residue acquision into one ptc_reference fucntion for mitochondria
# 2.use ptc_reference in consort with another structures's rrna to establish that structure's mapped residues
# 3.get the farthest pair of residues -- midpoint is the ptc in that structure

def trna_get_cterm_residues()->np.ndarray:
    rcsb_id, trna_id = REFERENCE_MITO_STRUCTURE_TRNA
    c:Chain = RibosomeOps(rcsb_id).biopython_structure()[0][trna_id]
    c_terminus:Residue = [*c][-1]
    return c_terminus.center_of_mass()

def mitorrna_ptc_residues(trna_cterm_pos:np.ndarray, mtrrna: Chain)->List[Residue]:
    atoms       = Selection.unfold_entities(mtrrna, "A")
    ns          = NeighborSearch(atoms)
    nbhd        = set()
    nearby_residues = ns.search(trna_cterm_pos, 10, "R")
    return nearby_residues


def get_ptc_mito(rcsb_id:str):
    trna_Cterm          = trna_get_cterm_residues()
    mtRRNA_src_aaid     = RibosomeOps(rcsb_id).get_LSU_rRNA().auth_asym_id
    mtRRNA_src_chain:Chain    = RibosomeOps(rcsb_id).biopython_structure()[0][mtRRNA_src_aaid]
    ress                = mitorrna_ptc_residues(trna_Cterm,mtRRNA_src)

    mtRRNA_target:Chain = RibosomeOps(target_rcsb_id).biopython_structure()[0][mtRRNA_target_aaid]
    _,_,motifs = map_motifs(BiopythonChain( mtRRNA_src ), BiopythonChain(mtRRNA_target), [ResidueSummary.from_biopython_residue(r) for r in ress], 'mt16SrRNA', True)
    x:Residue
    pprint(sorted(motifs, key=lambda x: x.get_id()))


# T is a landmark with method project_into, project_from, data D and flag `present`
# basically assgin to every node of the taxonomy tree the the landmark with the data where there is one
# "project" from extant nodes to the rest preferring proximal nodes as sources
class GlobalTaxonomy[T]():
    ...

def get_rcsb_ids():
    return rcsb_ids


trna_Cterm          = trna_get_cterm_residues()
mtRRNA_src:Chain    = RibosomeOps(src_rcsb_id).biopython_structure()[0][mtRRNA_src_aaid]
ress                = mitorrna_ptc_residues(trna_Cterm,mtRRNA_src)
mtRRNA_target:Chain = RibosomeOps(target_rcsb_id).biopython_structure()[0][mtRRNA_target_aaid]
_,_,motifs = map_motifs(BiopythonChain( mtRRNA_src ), BiopythonChain(mtRRNA_target), [ResidueSummary.from_biopython_residue(r) for r in ress], 'mt16SrRNA', True)
x:Residue
pprint(sorted(motifs, key=lambda x: x.get_id()))