import logging
import os
from rbxz_bend.db.ribosomexyz import ribosomexyzDB
from ribctl.lib.types import RibosomeStructure
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA, get_ribxz_logger
from ribctl.lib.mod_extract_bsites  import struct_ligand_ids, struct_liglike_ids, render_ligand, render_liglike_polymer
from ribctl.lib import utils

class App:
    def __init__(self):
        self.db_connection = db_connection
    
    def render_all_ligands(self):
        for file in os.listdir(RIBETL_DATA):
            if len(file) == 4 and os.path.isdir(os.path.join(RIBETL_DATA,file)) :
                PDBID = file
                logger = get_ribxz_logger('computations', __name__)

                _structure_cif_handle :Structure = utils.open_structure(PDBID,'cif')  # type: ignore
                struct_profile_handle       = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json'))  # type: ignore

                liglike_polys = struct_liglike_ids(struct_profile_handle)
                ligands       = struct_ligand_ids(PDBID, struct_profile_handle)

                print("Found ligands in {}: {}.".format(PDBID, ligands))
                for lig in ligands:
                    print(lig)
                print("Found ligandlike polymers in {}:".format(PDBID))
                for polyref in liglike_polys:
                    print(polyref)

                print("Rendering to file.")

                for polyref in liglike_polys:

                    try:
                        render_liglike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, True)
                        logger.info("Rendered liglike polymer {} in {}.".format(polyref.auth_asym_id, polyref.parent_rcsb_id))

                    except Exception as e:
                        logger.info("Error rendering liglike polymer {} in {}.".format(polyref.auth_asym_id, polyref.parent_rcsb_id))
                        logger.error(e)

                for l in ligands:

                    try:
                        render_ligand(PDBID, l[0], _structure_cif_handle, True)
                        logger.info("Rendered ligand {} in {}.".format(l[0], PDBID))

                    except Exception as e:
                        logger.info("Error rendering ligand {} in {}.".format(l[0], PDBID))
                        logger.error(e)

db_connection =  ribosomexyzDB(uri=NEO4J_URI,
                          password=NEO4J_PASSWORD,
                          user=NEO4J_USER)

ribosomexyzApp = App()

