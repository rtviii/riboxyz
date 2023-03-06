import json
import os
from neo4j import Driver
from pydantic import parse_obj_as
from ribctl.lib import RIBETL_DATA
from ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.struct_render_thumbnail import render_thumbnail
from ribctl.lib.struct_rcsb_api import process_pdb_record
from ribctl.lib.struct_split_rename import split_rename
from ribctl.lib.struct_extract_bsites import get_ligands, get_liglike_polymers, render_ligand, render_liglike_polymer





class RibosomeAssets():
    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()
# 
    def _envcheck(self):
        if not RIBETL_DATA:
            raise Exception(
                "RIBETL_DATA environment variable not set. Cannot access assets.")

    def _dir_path(self):
        self._envcheck()
        return f"{RIBETL_DATA}/{self.rcsb_id}"

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

    def biopython_structure(self):
        return open_structure(self.rcsb_id, 'cif')

    def chains_dir(self):
        self._envcheck()
        return f"{self._dir_path()}/CHAINS"

    def _png_thumbnail_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/_ray_{self.rcsb_id}.png"

    def save_json_profile(self, filepath: str, profile: dict):
        with open(filepath, "w") as f:
            json.dump(profile, f)

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def _verify_dir_exists(self):
        if not os.path.exists(self._dir_path()):
            os.umask(0)
            os.makedirs(self._dir_path(), 0o777)

    def _verify_cif(self, overwrite: bool = False) -> bool:
        if overwrite:
            download_unpack_place(self.rcsb_id)
            return True
        else:
            if os.path.exists(self._cif_filepath()):
                return True
            else:
                return False

    def _verify_cif_modified(self, overwrite: bool = False) -> bool:
        if overwrite:
            split_rename(self.rcsb_id)
            return True
        else:
            return os.path.exists(self._cif_modified_filepath())

    def _verify_json_profile(self, overwrite: bool = False) -> bool:
        if overwrite:
            ribosome = process_pdb_record(self.rcsb_id)
            if not parse_obj_as(RibosomeStructure, ribosome):
                raise Exception("Invalid ribosome structure profile.")

            self.save_json_profile(self._json_profile_filepath(), ribosome.dict())
            print(
                f"Saved structure profile:\t{self._json_profile_filepath()}")
            return True
        else:
            if os.path.exists(self._json_profile_filepath()):
                return True
            return False

    def _verify_png_thumbnail(self, overwrite: bool = False) -> bool:
        if overwrite:
            print("Obtaning thumbnail...")
            render_thumbnail(self.rcsb_id)
            return True
        else:
            if os.path.exists(self._png_thumbnail_filepath()):
                return True
            else:
                return False

    def _verify_chains_dir(self):
        split_rename(self.rcsb_id)

    def _verify_ligads_and_ligandlike_polys(self, overwrite: bool = False):

        def ligand_path(chem_id):            return os.path.join(self._dir_path(), f"LIGAND_{chem_id.upper()}.json")
        def liglike_poly_path(auth_asym_id): return os.path.join(self._dir_path(), f"POLYMER_{auth_asym_id.upper()}.json")

        ligands             = get_ligands(self.rcsb_id, self.json_profile())
        ligandlike_polymers = get_liglike_polymers(self.json_profile())

        _flag = True

        for ligand in ligands:
            if not os.path.exists(ligand_path(ligand[0])):
                _flag = False
                render_ligand(
                    self.rcsb_id, ligand[0], self.biopython_structure(), overwrite)

        for ligandlike_poly in ligandlike_polymers:
            if not os.path.exists(liglike_poly_path(ligandlike_poly.auth_asym_id)):
                _flag = False
                render_liglike_polymer(
                    self.rcsb_id, ligandlike_poly.auth_asym_id, self.biopython_structure(), overwrite)

        return _flag
