
import os
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB import PDBIO, Select
import numpy as np
from __archive.scripts.pymol_visualtion import extract_chains
from ribctl import EXIT_TUNNEL_WORK, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass
from pymol import cmd

RCSB_ID   = "4UG0"


# def extract_chains(rcsb_id: str):
#     outpath           = '{}/{}/{}_vestibule.mmcif'.format(EXIT_TUNNEL_WORK, RCSB_ID,RCSB_ID)
#     eukaryotic_chains = [ 'uL4', 'uL22', 'eL32', 'uL24', ]
#     auth_asym_ids     = []
#     for c  in eukaryotic_chains:
#         c_ = RibosomeAssets(rcsb_id).get_chain_by_polymer_class(c)
#         auth_asym_ids.append(c_.auth_asym_id)


#     cmd.load(path)
#     cmd.extract('crown', 'c. {}'.format('+'.join(auth_asym_ids)))
#     cmd.save(outpath, 'crown')
#     print("Wrote ", outpath)
#     return auth_asym_ids


def extract_chains_by_auth_asym_id(rcsb_id: str, chain_auth_asym_ids:list, outpath:str):
    path    = os.path.join(RIBETL_DATA, rcsb_id, "{}.cif".format(rcsb_id))
    cmd.load(path)
    cmd.extract('crown', 'c. {}'.format('+'.join(chain_auth_asym_ids)))
    cmd.save(outpath, 'crown')
    print("Wrote ", outpath)


# 