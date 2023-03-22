from ribctl.lib.mod_transpose_bsites import init_transpose_ligand, open_bsite
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib import utils
from ribctl.lib.mod_extract_bsites import save_ligand, save_ligandlike_polymer, save_ligandlike_polymer, struct_ligand_ids, struct_liglike_ids

PDBID = "3J7Z"

def repair_ligands (PDBID):
    _structure_cif_handle = utils.open_structure(PDBID,'cif')
    struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

    liglike_polys = struct_liglike_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(PDBID, struct_profile_handle)

    for polyref in liglike_polys:
        save_ligandlike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, True)

    for l in ligands:
        save_ligand(PDBID, l, _structure_cif_handle, True)


repair_ligands(PDBID)

# _structure_cif_handle = utils.open_structure(PDBID,'cif')
# struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

# liglike_polys = struct_liglike_ids(struct_profile_handle)
# ligands       = struct_ligand_ids(PDBID, struct_profile_handle)




# bsite          = open_bsite(utils.ligand_path(PDBID, 'ERY', 'LIGAND'))
# target         = '5afi'
# target_profile = RibosomeStructure.parse_obj(utils.open_structure(target,'json')  )
# prediction = init_transpose_ligand(target_profile, bsite)