import json
import os
from pprint import pprint
from pydantic import BaseModel, parse_obj_as
import gzip
import os
import requests
from utils import download_unpack_place, open_structure
from types_ribosome import RibosomeStructure
from process_structure import process_pdb_record
from render_thumbnail import render_thumbnail
from split_rename import split_rename
from extract_bsites import get_ligands, get_liglike_polymers, render_ligand, render_liglike_polymer
from Bio.PDB.Structure import Structure
from Bio.PDB import FastMMCIFParser
import typing

RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))

class RibosomeAssets():
    rcsb_id: str
    def __init__(self, rcsb_id:str) -> None:
        self.rcsb_id= rcsb_id

    def _envcheck(self):
        if not os.environ.get("RIBETL_DATA"):
            raise Exception(
                "RIBETL_DATA environment variable not set. Cannot access assets.")

    def _dir_path(self):
        self._envcheck()
        return f"{os.environ.get('RIBETL_DATA')}/{self.rcsb_id}"

    def _cif_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.cif"

    def _cif_modified_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}_modified.cif"

    def _json_profile_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.json"

    def json_profile(self):
        with open(self._json_profile_filepath(), "r") as f:
            return json.load(f)

    def biopython_sturcture(self):
        return open_structure(self.rcsb_id, 'cif')

    def chains_dir(self):
        self._envcheck()
        return f"{self._dir_path()}/CHAINS"

    def _png_thumbnail_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/_ray_{self.rcsb_id}.png"
      
    def save_json_profile(self,filepath:str, profile:dict):
        with open(filepath, "w") as f:
            json.dump(profile, f)

    def _verify_cif(self, obtain: bool = False) -> bool:
        if os.path.exists(self._cif_filepath()):
            return True
        else:
            if obtain:
                download_unpack_place(self.rcsb_id)
                return True
            else:
                return False

    def _verify_cif_modified(self, obtain: bool = False) -> bool:
        if os.path.exists(self._cif_modified_filepath()):
            return True
        else:
            if obtain:
                split_rename(self.rcsb_id)
                return True
            else:
                return False

    def _verify_json_profile(self, obtain: bool = False) -> bool:
        if os.path.exists(self._json_profile_filepath()):
            return True
        else:
            if obtain:
                ribosome = process_pdb_record(self.rcsb_id)
                if not parse_obj_as(RibosomeStructure, ribosome):
                    raise Exception("Invalid ribosome structure profile.")
                  
                self.save_json_profile(self._json_profile_filepath(), ribosome)
                print(f"Saved structure profile:\t{self._json_profile_filepath()}")
                return True
            else:
                return False

    def _verify_png_thumbnail(self, obtain: bool = False) -> bool:
        if os.path.exists(self._png_thumbnail_filepath()):
            return True
        else:
            if obtain:
                print("Obtaning thumbnail...")
                render_thumbnail(self.rcsb_id)

                return True
            else:
                return False
    
    def _verify_chains_dir(self, obtain: bool = False):
        split_rename(self.rcsb_id)
        
    def _verify_ligads_and_ligandlike_polys(self, obtain:bool=False):

        ligand_path         = lambda chem_id: os.path.join(self._dir_path(), f"LIGAND_{chem_id.upper()}.json")
        liglike_poly_path   = lambda auth_asym_id: os.path.join(self._dir_path(), f"POLYMER_{auth_asym_id.upper()}.json")

        ligands             = get_ligands(self.rcsb_id, self.json_profile())
        ligandlike_polymers = get_liglike_polymers(self.json_profile())

        _flag               = True

        for ligand in ligands:
            if not os.path.exists(ligand_path(ligand[0])):
                _flag = False
                render_ligand(self.rcsb_id, ligand[0],self.biopython_sturcture(), obtain)

        for ligandlike_poly in ligandlike_polymers:
            if not os.path.exists(liglike_poly_path(ligandlike_poly.auth_asym_id)):
                _flag = False
                render_liglike_polymer(self.rcsb_id, ligandlike_poly.auth_asym_id,self.biopython_sturcture(), obtain)

        return _flag

new = RibosomeAssets("3J7Z")

new._verify_cif(True)
new._verify_json_profile(True)
new._verify_cif_modified(True)
new._verify_ligads_and_ligandlike_polys(True)
new._verify_chains_dir(True)

      
        