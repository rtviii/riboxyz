from enum import   auto
import json
import os
from pprint import pprint
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from loguru import logger
import requests
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass,  RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA

class StructureAssetPaths:
    rcsb_id:str
    def __init__(self, rcsb_id) -> None:
        self.rcsb_id = rcsb_id
        pass

    def binding_site(self, chemId:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}.json"

    def binding_site_prediction(self, chemId:str, source_struct:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}_PREDICTION_VIA_{source_struct.upper()}.json"
    
    @property
    def dir(self):
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    @property
    def cif(self):
        return f"{self.dir}/{self.rcsb_id}.cif"

    @property
    def ptc(self):
        return os.path.join(self.dir, "{}_PTC.json".format(self.rcsb_id) )

    @property
    def profile(self):
        return os.path.join(self.dir, f"{self.rcsb_id}.json")

    @property
    def chains_dir(self):
        return f"{self.dir}/CHAINS"

    @property
    def classification_report(self):
        return os.path.join(self.dir, f"classification_report_{self.rcsb_id}.json")

    @property
    def thumbnail(self):
        return f"{self.dir}/{self.rcsb_id}.png"

class StructureAssets:

    rcsb_id: str
    paths  : StructureAssetPaths

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id
        self.paths   = StructureAssetPaths(rcsb_id)
    
    def biopython_structure(self)-> Structure:
        cifpath = StructureAssetPaths(self.rcsb_id).cif
        return FastMMCIFParser(QUIET=True).get_structure(self.rcsb_id, cifpath)

    def profile(self) -> RibosomeStructure:
        with open(self.paths.profile, "r") as f:
            return RibosomeStructure.model_validate(json.load(f))

    def ptc(self) -> PTCInfo:
        with open(self.paths.ptc, "r") as infile:
            _ = json.load(infile)
            return PTCInfo.model_validate(_)

    def biopython_get_chain(self, auth_asym_id: str) -> Chain:
        return self.biopython_structure().child_dict[0].child_dict[auth_asym_id]

    def acquire_all_assets(self):
        ...


    def _verify_dir_exists(self):
        if not os.path.exists(self.paths.dir):
            os.umask(0)
            os.makedirs(self.paths.dir, 0o755)