from pprint import pprint
from typing import List, Tuple, TypeVar
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
import numpy as np
from ribctl.etl.etl_assets_ops import RibosomeOps
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection

from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import SequenceMappingContainer
from ribctl.lib.schema.types_binding_site import ResidueSummary
from scipy.spatial.distance import pdist, squareform

mitochondrial_rcsb_ids = [
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
    "8QSJ",
]



def find_closest_pair(points:np.ndarray):
    points = np.asarray(points)
    if len(points) < 2:
        raise ValueError("Array must contain at least 2 points")
    
    distances = pdist(points)
    distance_matrix = squareform(distances)
    i, j = np.triu_indices(len(points), k=1)
    min_idx = np.argmin(distances)
    point1_idx = i[min_idx]
    point2_idx = j[min_idx]
    
    closest_point1 = points[point1_idx]
    closest_point2 = points[point2_idx]
    min_distance   = distances[min_idx]
    
    return closest_point1, closest_point2, min_distance
# TODO:
# 1.combine the cterm/residue acquision into one ptc_reference fucntion for mitochondria
# 2.use ptc_reference in consort with another structures's rrna to establish that structure's mapped residues
# 3.get the farthest pair of residues -- midpoint is the ptc in that structure
# 4.feed to mesh acq

REFERENCE_MITO_STRUCTURE_TRNA_RRNA = ("7A5F", "24", "A3")
ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = REFERENCE_MITO_STRUCTURE_TRNA_RRNA


def PTC_reference_mito() -> List[Residue]:
    mmcif_struct = RibosomeOps(ref_rcsb_id).biopython_structure()[0]

    def trna_cterm_pos() -> np.ndarray: 
        trnaChain                     : Chain = mmcif_struct[ref_trna_aaid]
        c_terminus                    : Residue = [*trnaChain][-1]
        return c_terminus.center_of_mass()

    mtrrna          = mmcif_struct[ref_rrna_aaid]
    atoms           = Selection.unfold_entities(mtrrna, "A")
    ns              = NeighborSearch(atoms)
    nearby_residues = ns.search(trna_cterm_pos(), 10, "R")

    return nearby_residues


def get_ptc_mito(target_rcsb_id: str)->Tuple[np.ndarray ,list[Residue]]:
    """
    Get PTC in @target_rcsb_id by way of mapping a reference PTC in a given mitochondrial structure
    """
    mmcif_struct_src = RibosomeOps(ref_rcsb_id).biopython_structure()[0]
    mmcif_struct_tgt = RibosomeOps(target_rcsb_id).biopython_structure()[0]

    mtRRNA_tgt_aaid = RibosomeOps(target_rcsb_id).get_LSU_rRNA().auth_asym_id

    rrna_src = mmcif_struct_src[ref_rrna_aaid]
    rrna_tgt = mmcif_struct_tgt[mtRRNA_tgt_aaid]

    ref_residues = PTC_reference_mito()

    _, _, motifs = map_motifs(
        SequenceMappingContainer(rrna_src),
        SequenceMappingContainer(rrna_tgt),
        [ResidueSummary.from_biopython_residue(r) for r in ref_residues],
        "mt16SrRNA",
        True )

    (p1,p2,dist) = find_closest_pair([r.center_of_mass()  for r in motifs])
    return (p1+p2)/2, motifs


# T is a landmark with method project_into, project_from, data D and flag `present`
# basically assgin to every node of the taxonomy tree the the landmark with the data where there is one
# "project" from extant nodes to the rest preferring proximal nodes as sources
# class GlobalTaxonomy[T](): ...
