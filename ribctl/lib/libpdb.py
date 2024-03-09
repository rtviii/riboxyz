
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO, Select
import numpy as np
from __archive.scripts.pymol_visualtion import extract_chains
from ribctl import EXIT_TUNNEL_WORK, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass
from pymol import cmd



def extract_chains_by_auth_asym_id(rcsb_id: str, chain_auth_asym_ids:list, outpath:str):
    path    = os.path.join(RIBETL_DATA, rcsb_id, "{}.cif".format(rcsb_id))
    cmd.load(path)
    cmd.extract('crown', 'c. {}'.format('+'.join(chain_auth_asym_ids)))
    cmd.save(outpath, 'crown')
    print("Wrote ", outpath)
