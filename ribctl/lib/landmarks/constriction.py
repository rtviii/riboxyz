from Bio.PDB.Chain import Chain
import numpy as np
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.landmarks.ptc_via_doris import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import find_closest_pair_two_sets, midpoint
from ribctl.lib.utils import download_unpack_place
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass, PolynucleotideClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure, RibosomeStructureMetadata, )
from ribctl.logs.loggers import get_etl_logger

def get_constriction(rcsb_id: str)->np.ndarray:
    ro               = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile.mitochondrial
    if is_mitochondrial:
        uL4  = ro.get_poly_by_polyclass('uL4m')
        uL22 = ro.get_poly_by_polyclass('uL22m')
    else:
        uL4  = ro.get_poly_by_polyclass('uL4')
        uL22 = ro.get_poly_by_polyclass('uL22')

    if uL4 is None or uL22 is None:
        raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))

    structure = ro.assets.biopython_structure()

    uL4_c :Chain = structure[0][uL4.auth_asym_id]
    uL22_c :Chain = structure[0][uL22.auth_asym_id]

    uL4_coords  = [(r.center_of_mass() ) for r in uL4_c.child_list]
    uL22_coords = [(r_.center_of_mass() ) for r_ in uL22_c.child_list]

    return midpoint(*find_closest_pair_two_sets(uL4_coords, uL22_coords))