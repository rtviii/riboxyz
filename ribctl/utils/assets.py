import json
import os
from pydantic import BaseModel, parse_obj_as
import gzip
import os
import requests
from render_thumbnail import render_thumbnail
from ribctl.types_ribosome import RibosomeStructure
from ribctl.utils.process_structure import process_pdb_record
from split_rename import split_rename
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser
import typing
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))


def download_unpack_place(struct_id: str) -> None:

    BASE_URL = "http://files.rcsb.org/download/"
    FORMAT = ".cif.gz"

    structid = struct_id.upper()
    url = BASE_URL + structid + FORMAT
    compressed = requests.get(url).content
    decompressed = gzip.decompress(compressed)

    destination_chains = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        "CHAINS")

    if not os.path.exists(destination_chains):
        os.mkdir(destination_chains)
        print(f"Created directory {destination_chains}.")

    structfile = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        structid + ".cif")

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


class RibosomeAssets(BaseModel):
    rcsb_id: str

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id

    def envcheck(self):
        if not os.environ.get("RIBETL_DATA"):
            raise Exception(
                "RIBETL_DATA environment variable not set. Cannot access assets.")

    def dir_path(self):
        self.envcheck()
        return f"{os.environ.get('RIBETL_DATA')}/{self.rcsb_id}"

    def cif_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}.cif"

    def cif_modified_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}_modified.cif"

    def json_profile_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/{self.rcsb_id}.json"

    def chains_folder(self):
        self.envcheck()
        return f"{self.dir_path()}/CHAINS"

    def png_thumbnail_filepath(self):
        self.envcheck()
        return f"{self.dir_path()}/_ray_{self.rcsb_id}.png"
      
    @staticmethod
    def save_json_profile(filepath:str, profile:dict):
        with open(filepath, "w") as f:
            json.dump(profile, f)

    def __verify_cif(self, obtain: bool = False) -> bool:
        if os.path.exists(self.cif_filepath()):
            return True
        else:
            if obtain:
                download_unpack_place(self.rcsb_id)
                return True
            else:
                return False

    def __verify_cif_modified(self, obtain: bool = False) -> bool:
        if os.path.exists(self.cif_modified_filepath()):
            return True
        else:
            if obtain:
                split_rename(self.rcsb_id)
                return True
            else:
                return False

    def __verify_json_profile(self, obtain: bool = False) -> bool:
        if os.path.exists(self.json_profile_filepath()):
            return True
        else:
            if obtain:
                ribosome = process_pdb_record(self.rcsb_id)
                if not parse_obj_as(RibosomeStructure, ribosome):
                    raise Exception("Invalid ribosome structure profile.")
                  
                self.save_json_profile(self.json_profile_filepath(), ribosome)
                print(f"Saved structure profile:\t{self.json_profile_filepath()}")
                return True
            else:
                return False

    def __verify_png_thumbnail(self, obtain: bool = False) -> bool:
        if os.path.exists(self.png_thumbnail_filepath()):
            return True
        else:
            if obtain:
                print("Obtaning thumbnail...")
                render_thumbnail(self.rcsb_id)

                return True
            else:
                return False
    
    def __verify_chains_folder(self, obtain: bool = False) -> bool:
        if os.path.exists(self.chains_folder()):
            return True
        else:
            if obtain:
                self.__verify_cif_modified(True)
                return True
            else:
                return False

    def verify_ligads_and_ligandlike_polys(self):

        ...
      
        