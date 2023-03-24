import os
from pprint import pprint
from ribctl.lib.types.ligands.types_binding_site import BindingSite, LigandPrediction
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand, open_bsite
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib import RIBETL_DATA, utils
from ribctl.lib.mod_extract_bsites import save_ligand, save_ligandlike_polymer, save_ligandlike_polymer, struct_ligand_ids, struct_liglike_ids
from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_polyclass

PDBID1     = "7K00"
PDBID2     = "5AFI"
resrange   = (0,50)
poly_class = "5SrRNA"

def extract_bsites (rcsb_id):
    _structure_cif_handle = utils.open_structure(rcsb_id,'cif')
    struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(rcsb_id,'json')  )

    liglike_polys = struct_liglike_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(rcsb_id, struct_profile_handle)

    for polyref in liglike_polys:
        bsite_p:BindingSite = save_ligandlike_polymer(polyref.auth_asym_id, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'POLYMER_{polyref.auth_asym_id}.json')

        bsite_p.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

    for chemid in ligands:
        bsite_l = save_ligand( chemid, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'LIGAND_{chemid}.json')

        bsite_l.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

def predict_bsite(PDBID_target, ligand_id, ligand_type)->LigandPrediction:

    bsite          = open_bsite(utils.ligand_path(PDBID_target, ligand_id, ligand_type))
    target_profile = RibosomeStructure.parse_obj(utils.open_structure(PDBID_target,'json')  )
    prediction     = init_transpose_ligand(target_profile, bsite)

    return prediction

# (src_auth_asym_id, src_path, src_range,
#  tgt_auth_asym_id, tgt_path, tgt_range) = ranged_align_by_polyclass(PDBID1,PDBID2,resrange,poly_class)

# cif_str = pymol_super(
#     PDBID1,
#     src_range,
#     src_auth_asym_id,

#     PDBID2,
#     tgt_range,
#     tgt_auth_asym_id,
# )




# extract_bsites(PDBID)
# print(cif_str)
