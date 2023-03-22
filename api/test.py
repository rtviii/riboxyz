from pprint import pprint
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand, open_bsite
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib import utils
from ribctl.lib.mod_extract_bsites import save_ligand, save_ligandlike_polymer, save_ligandlike_polymer, struct_ligand_ids, struct_liglike_ids

# PDBID = "3J7Z"
PDBID = "5AFI"

def repair_ligands (PDBID):
    _structure_cif_handle = utils.open_structure(PDBID,'cif')
    struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

    liglike_polys = struct_liglike_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(PDBID, struct_profile_handle)

    for polyref in liglike_polys:
        save_ligandlike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, True)

    for l in ligands:
        save_ligand(PDBID, l, _structure_cif_handle, True)


# repair_ligands(PDBID)
# 
# _structure_cif_handle = utils.open_structure(PDBID,'cif')
# struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

# liglike_polys = struct_liglike_ids(struct_profile_handle)
# ligands       = struct_ligand_ids(PDBID, struct_profile_handle)




# bsite          = open_bsite(utils.ligand_path(PDBID, 'ERY', 'LIGAND'))
# bsite          = open_bsite(utils.ligand_path(PDBID, '7', 'POLYMER'))
# ---------
# drwxrwxrwx    3 root root 4.0K Feb  7 14:26 .
# drwxrwxrwx 1442 root root  48K Mar 19 00:38 ..
# -rwxrwxrwx    1 root root  26M Mar  6 12:28 5AFI.cif
# -rwxrwxrwx    1 root root  66K Mar 20 21:28 5AFI.json
# -rwxrwxrwx    1 root root  22M Oct 21  2021 5AFI_modified.cif
# drwxrwxrwx    2 root root 4.0K Oct 21  2021 CHAINS
# -rwxrwxrwx    1 root root 8.0K Feb  6 00:12 LIGAND_GDP.json
# -rwxrwxrwx    1 root root 6.2K Feb  6 00:12 LIGAND_KIR.json
# -rwxrwxrwx    1 root root  444 Feb  7 14:26 l-s-whole.md
# -rwxrwxrwx    1 root root  18K Feb  6 00:12 POLYMER_v.json
# -rwxrwxrwx    1 root root  19K Feb  6 00:12 POLYMER_w.json
# -rwxrwxrwx    1 root root  14K Feb  6 00:12 POLYMER_x.json
# -rwxrwxrwx    1 root root  22K Feb  6 00:12 POLYMER_y.json
# -rwxrwxrwx    1 root root  13K Feb  6 00:12 POLYMER_z.json
# -rwxrwxrwx    1 root root 421K Feb  2 19:45 _ray_5AFI.png


# bsite          = open_bsite(utils.ligand_path(PDBID, 'GDP', 'LIGAND'))
bsite          = open_bsite(utils.ligand_path(PDBID, 'w', 'POLYMER'))
target         = '7k00'
target_profile = RibosomeStructure.parse_obj(utils.open_structure(target,'json')  )
prediction = init_transpose_ligand(target_profile, bsite)

pprint("Prediction")
pprint(prediction)