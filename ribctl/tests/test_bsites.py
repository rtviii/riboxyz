import os
from pprint import pprint
from ribctl.lib.ribosome_types.types_binding_site import BindingSite, LigandPrediction
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand, open_bsite
from ribctl.lib.ribosome_types.types_ribosome import RibosomeStructure
from ribctl.lib import utils
from ribctl.lib.mod_extract_bsites import bsite_nonpolymeric_ligand, bsite_polymeric_factor, bsite_polymeric_factor, struct_ligand_ids, struct_polymeric_factor_ids
from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_polyclass

PDBID1     = "7K00"
PDBID2     = "5AFI"
resrange   = (0,50)
poly_class = "5SrRNA"

def extract_bsites (rcsb_id):
    _structure_cif_handle = utils.open_structure(rcsb_id,'cif')
    struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(rcsb_id,'json')  )

    liglike_polys = struct_polymeric_factor_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(rcsb_id, struct_profile_handle)

    for polyref in liglike_polys:
        bsite_p:BindingSite = bsite_polymeric_factor(polyref.auth_asym_id, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'POLYMER_{polyref.auth_asym_id}.json')

        bsite_p.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

    for chemid in ligands:
        bsite_l = bsite_nonpolymeric_ligand( chemid, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'LIGAND_{chemid}.json')

        bsite_l.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

def predict_bsite(PDBID_target, ligand_id, ligand_type)->LigandPrediction:

    bsite          = open_bsite(utils.ligand_path(PDBID_target, ligand_id, ligand_type))
    target_profile = RibosomeStructure.parse_obj(utils.open_structure(PDBID_target,'json')  )
    prediction     = init_transpose_ligand(target_profile, bsite)

    return prediction
