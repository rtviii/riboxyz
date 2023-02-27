import gzip
import json
import os
import typing
import requests
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser

RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))

def download_unpack_place(struct_id: str) -> None:

    BASE_URL = "http://files.rcsb.org/download/"
    FORMAT   = ".cif.gz"

    structid           = struct_id.upper()
    url                = BASE_URL + structid + FORMAT
    compressed         = requests.get(url).content
    decompressed       = gzip.decompress(compressed)

    destination_chains = os.path.join(
            os.environ["RIBETL_DATA"],
            structid,
            "CHAINS"
        )

    if not os.path.exists(destination_chains):
        os.mkdir(destination_chains)
        print(f"Created directory {destination_chains}.")

    structfile = os.path.join(
            os.environ["RIBETL_DATA"],
            structid,
            structid + ".cif"
        )

    with open(structfile, "wb") as f:
        f.write(decompressed)

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

    elif path_type == 'modified':
        with open(struct_path(pdbid, 'modified'), 'rb') as _:

            try:
                return FastMMCIFParser(QUIET=True).get_structure(pdbid, _)

            except Exception as e:
                return f"\033[93m Parser Error in structure {pdbid} \033[0m : {e}"