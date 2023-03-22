import os
from pprint import pprint
from ribctl.lib.types.ligands.types_binding_site import BindingSite, LigandPrediction
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand, open_bsite
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib import RIBETL_DATA, utils
from ribctl.lib.mod_extract_bsites import save_ligand, save_ligandlike_polymer, save_ligandlike_polymer, struct_ligand_ids, struct_liglike_ids
from ribctl.lib.mod_superimpose import ranged_super_by_polyclass

# PDBID = "3J7Z"
PDBID = "5AFI"

def extract_bsites (PDBID):
    _structure_cif_handle = utils.open_structure(PDBID,'cif')
    struct_profile_handle = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json')  )

    liglike_polys = struct_liglike_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(PDBID, struct_profile_handle)

    for polyref in liglike_polys:
        bsite_p:BindingSite = save_ligandlike_polymer(polyref.auth_asym_id, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, PDBID.upper(), f'POLYMER_{polyref.auth_asym_id}.json')
        bsite_p.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)


    for chemid in ligands:
        bsite_l = save_ligand( chemid, _structure_cif_handle)
        outfile_json = os.path.join(RIBETL_DATA, PDBID.upper(), f'LIGAND_{chemid}.json')
        bsite_l.save(outfile_json)
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

def predict_bsite(PDBID_target, ligand_id, ligand_type)->LigandPrediction:
    bsite          = open_bsite(utils.ligand_path(PDBID_target, ligand_id, ligand_type))
    target_profile = RibosomeStructure.parse_obj(utils.open_structure(PDBID_target,'json')  )
    prediction     = init_transpose_ligand(target_profile, bsite)

    return prediction

print(ranged_super_by_polyclass('3j7z','5afi',(0,20),'uL4'))

# extract_bsites(PDBID)
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
