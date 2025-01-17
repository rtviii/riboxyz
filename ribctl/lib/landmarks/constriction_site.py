from Bio.PDB.Chain import Chain
import numpy as np
from ribctl.lib.types.polymer.base import CytosolicProteinClass, MitochondrialProteinClass
from ribctl.lib.utils import find_closest_pair_two_sets, midpoint
from ribctl.ribosome_ops import RibosomeOps

def get_constriction(rcsb_id: str)->np.ndarray:
    ro               = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile.mitochondrial
    if is_mitochondrial:
        uL4  = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL4m)
        uL22 = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL22m)
    else:
        uL4  = ro.get_poly_by_polyclass(CytosolicProteinClass.uL4)
        uL22 = ro.get_poly_by_polyclass(CytosolicProteinClass.uL22)
    if uL4 is None or uL22 is None:
        raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))

    structure = ro.assets.biopython_structure()
    
    uL4_c  :Chain = structure[0][uL4.auth_asym_id]
    uL22_c :Chain = structure[0][uL22.auth_asym_id]

    uL4_coords  = [(r.center_of_mass() ) for r in uL4_c.child_list]
    uL22_coords = [(r_.center_of_mass() ) for r_ in uL22_c.child_list]

    return midpoint(*find_closest_pair_two_sets(uL4_coords, uL22_coords))
