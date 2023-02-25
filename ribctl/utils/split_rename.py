import argparse
import json
import os
import typing
from Bio.PDB.Structure import Structure
import itertools
import itertools
import numpy as np
from Bio.PDB import MMCIF2Dict, MMCIFIO, FastMMCIFParser
import gemmi
flatten = itertools.chain.from_iterable
n1 = np.array
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))


def __get_dict(path: str,) -> dict:
    return MMCIF2Dict.MMCIF2Dict(path)


def make_nom_dict(profile):
    nomdict = {}
    for i in [*profile['rnas'], *profile['proteins']]:
        nomdict[i['auth_asym_id']] = i['nomenclature']
    return nomdict


def struct_path(pdbid: str, pftype: typing.Literal["cif", "json", "modified"]):
    if pftype == 'cif':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}.cif")
    elif pftype == 'json':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}.json")
    elif pftype == 'modified':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}_modified.cif")
    else:
        raise ValueError("Invalid path type. Must be 'cif', 'json', or 'modified' ")


def open_structure(pdbid: str, path_type: typing.Literal[ "cif", "json", "modified"])->Structure|typing.Any:
    pdbid = pdbid.upper()
    if path_type == 'cif':
        cifpath = struct_path(pdbid, 'cif')
        try:
            return FastMMCIFParser(QUIET=True).get_structure(pdbid, cifpath)
        except Exception as e:
            return f"\033[93m Parser Error in structure {pdbid} \033[0m : {e}"

    if path_type == 'json':
        with open(struct_path(pdbid, 'json'), 'rb') as _:
            return json.load(_)


def __inject_dict(pdbid: str):

    cifpath = os.path.join(RIBETL_DATA, pdbid, pdbid+".cif")
    cifmodified = os.path.join(RIBETL_DATA, pdbid, pdbid+"_modified.cif")
    if os.path.exists(cifmodified):
        print("Skipping modified struct. Exists: \t", cifmodified)
        return
    structprofile = open_structure(pdbid, 'json')

    doc = gemmi.cif.read_file(cifpath)
    block = doc.sole_block()
    loop = block.init_loop('_ribosome_nomenclature.', [
                           'entity_poly.pdbx_strand_id', 'polymer_class'])

    nomd = make_nom_dict(structprofile)
    for i in make_nom_dict(structprofile).items():
        loop.add_row([i[0], gemmi.cif.quote(
            "unclassified" if len(i[1]) == 0 else i[1][0])])

    doc.write_file(cifmodified)
    print("\033[91mWrote {} \033[0m".format(cifmodified))


def __process_chains(pdbid: str):
    io = MMCIFIO()
    if not os.path.exists(os.path.join(RIBETL_DATA, pdbid, 'CHAINS')):
        os.mkdir(os.path.join(RIBETL_DATA, pdbid, 'CHAINS'))

    structprofile = open_structure(pdbid, 'json')
    struct_cif: Structure = open_structure(pdbid, 'cif')

    model = gemmi.read_structure(os.path.join(
        RIBETL_DATA, pdbid, f'{pdbid}.cif'))[0]
    nomd = make_nom_dict(structprofile)

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
                print("Saved ", destination)
            except:
                ...
            io.save(destination)
        else:
            print("Skipping chain. Exists: \t", destination)


def split_rename(pdbid: str):
    pdbid = pdbid.upper()
    __inject_dict(pdbid)
    __process_chains(pdbid)


if __name__ == '__main__':
    parser = argparse. ArgumentParser(
        description='Split structure into constituent polymers and inject new nomencalture into the .cif file')
    parser. add_argument("-s", "--structure", type=str,
                         help="RCSB ID of structure to process")

    args = parser.parse_args()
    pdbid = args.structure

    if not pdbid:
        raise ValueError("Provide structure ID with -s arg")

    else:
        __inject_dict(pdbid.upper())
        __process_chains(pdbid.upper())
