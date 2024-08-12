from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
import gzip
import os
import requests

async def download_unpack_place(struct_id: str) -> None:

    BASE_URL = "http://files.rcsb.org/download/"
    FORMAT   = ".cif.gz"

    structid     = struct_id.upper()
    url          = BASE_URL + structid + FORMAT
    compressed   = requests.get(url).content
    decompressed = gzip.decompress(compressed)

    structfile = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        structid + ".cif")

    with open(structfile, "wb") as f:
        f.write(decompressed)


def open_structure(pdbid: str) -> Structure :
    ...
