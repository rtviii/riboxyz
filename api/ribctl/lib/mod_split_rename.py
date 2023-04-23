import argparse
import os
import itertools
import itertools
import numpy as np
import gemmi
from Bio.PDB import MMCIF2Dict, MMCIFIO
from Bio.PDB.Structure import Structure
from ribctl.lib.utils import open_structure

flatten     = itertools.chain.from_iterable
n1          = np.array
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))

def __get_dict(path: str,) -> dict:
    return MMCIF2Dict.MMCIF2Dict(path)

def __make_nom_dict(profile):
    nomdict = {}
    for i in [*profile['rnas'], *profile['proteins']]:
        nomdict[i['auth_asym_id']] = i['nomenclature']
    return nomdict

def __inject_dict(pdbid: str):

    cifpath     = os.path.join(RIBETL_DATA, pdbid, pdbid+".cif")
    cifmodified = os.path.join(RIBETL_DATA, pdbid, pdbid+"_modified.cif")

    if os.path.exists(cifmodified):
        print("Skipping modified struct. Exists: \t", cifmodified)
        return

    structprofile = open_structure(pdbid, 'json')

    doc = gemmi.cif.read_file(cifpath)
    block = doc.sole_block()
    loop = block.init_loop('_ribosome_nomenclature.', [
                           'entity_poly.pdbx_strand_id', 'polymer_class'])

    nomd = __make_nom_dict(structprofile)
    for i in __make_nom_dict(structprofile).items():
        loop.add_row([i[0], gemmi.cif.quote(
            "unclassified" if len(i[1]) == 0 else i[1][0])])

    doc.write_file(cifmodified)
    print("\033[91mWrote {} \033[0m".format(cifmodified))

def __process_chains(pdbid: str):
    io = MMCIFIO()
    if not os.path.exists(os.path.join(RIBETL_DATA, pdbid, 'CHAINS')):
        os.mkdir(os.path.join(RIBETL_DATA, pdbid, 'CHAINS'))

    structprofile = open_structure(pdbid, 'json')
    struct_cif: Structure    = open_structure(pdbid, 'cif')

    model = gemmi.read_structure(os.path.join(
        RIBETL_DATA, pdbid, f'{pdbid}.cif'))[0]
    nomd = __make_nom_dict(structprofile)

    for chain in struct_cif[0].child_list:
        destination = os.path.join(
            RIBETL_DATA, pdbid, 'CHAINS', '{}_STRAND_{}.cif'.format(pdbid, chain.get_id()))
        if not os.path.exists(destination):
            io.set_structure(chain)
            io.save(destination)
            cdict = __get_dict(destination)

            try:
                if len(nomd[chain.get_id()]) < 1:
                    cdict['data_'] = f'{pdbid}_{chain.get_id()}'
                else:
                    cdict['data_'] = f'{pdbid}_{chain.get_id()}_{nomd[chain.get_id()][0]}'
                io.set_dict(cdict)
                print("Saved: ", destination)
            except:
                ...
            io.save(destination)
        else:
            ...
            # print("Exists: ", destination)

async def split_rename(pdbid: str):
    pdbid = pdbid.upper()
    __inject_dict(pdbid)
    __process_chains(pdbid)


if __name__ == '__main__':
    parser = argparse. ArgumentParser(
        description='Split structure into constituent polymers and inject new nomencalture into the .cif file')
    parser.add_argument("-s", "--structure", type=str,
                         help="RCSB ID of structure to process")

    args  = parser.parse_args()
    pdbid = args.structure

    if not pdbid:
        raise ValueError("Provide structure ID with -s arg")
    else:
        __inject_dict(pdbid.upper())
        __process_chains(pdbid.upper())
