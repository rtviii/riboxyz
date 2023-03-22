from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib import utils
from ribctl.lib.struct_extract_bsites import render_liglike_polymer, render_ligand, struct_ligand_ids, struct_liglike_ids

PDBID = "3J7Z"

_structure_cif_handle = utils.open_structure(PDBID,'cif')
struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

liglike_polys = struct_liglike_ids(struct_profile_handle)
ligands       = struct_ligand_ids(PDBID, struct_profile_handle)

for polyref in liglike_polys:
    render_liglike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, True)

for l in ligands:
    render_ligand(PDBID, l[0], _structure_cif_handle, True)
