import os
from db.ribosomexyz import ribosomexyzDB
from ribctl.lib.ribosome_types import RibosomeStructure
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA
from ribctl.lib.mod_extract_bsites  import struct_ligand_ids, bsite_extrarbx_polymer, bsite_extrarbx_polymer
from ribctl.lib import utils

class App:
    def __init__(self):
        self.db_connection = db_connection
    

db_connection =  ribosomexyzDB(uri=NEO4J_URI,
                          password=NEO4J_PASSWORD,
                          user=NEO4J_USER)

ribosomexyzApp = App()

