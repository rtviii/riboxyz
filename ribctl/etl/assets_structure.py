import asyncio
from enum import   auto
import enum
import json
import os
from pprint import pprint
import subprocess
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from loguru import logger
import requests
from ribctl import AMINO_ACIDS_3_TO_1_CODE, ASSETS_PATH, CHAINSPLITTER_PATH, CLASSIFICATION_REPORTS
from ribctl.lib.libtax import PhylogenyNode, PhylogenyRank, Taxid
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from ribctl.lib.landmarks.ptc_via_doris import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.utils import download_unpack_place
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass, PolynucleotideClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA
from ribctl.logs.loggers import get_etl_logger


"""
Interface to the structural assets corresponding to a single entry in the database and a single subdir of `RIBETL_DATA`.
Example layout for a given structure:
    RIBETL_DATA/
        ├── 4UG0.cif
        ├── 4UG0.json
        ├── 4UG0_PTC.json
        ├── TUNNELS
        │   ├── 4UG0_convex_hull.npy
        │   ├── 4UG0_normal_estimated_surf.ply
        │   ├── 4UG0_poisson_recon.ply
        │   ├── 4UG0_poisson_recon_ascii.ply
        │   ├── 4UG0_spheres_expanded_pointset.npy
        │   ├── 4UG0_tunnel_atoms_bbox.json
        │   └── tunnel_4UG0.csv
        └── classification_report_4UG0.json
"""


class AssetClass(enum.StrEnum):
    profile   = auto()
    mmcif     = auto()

    thumbnail = auto()

    # landmarks
    ptc         = auto()
    npet        = auto()
    rna_helices = auto()
    trna_sites  = auto()

    @staticmethod
    def from_str(_:str):
        assets = [status for status in dir( AssetClass) if not status.startswith('_')]
        if _ in assets:
            return getattr(AssetClass, _)
        return None

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

    @property
    def tunnel_dir(self):
        if not os.path.exists(os.path.join(self.dir, "TUNNELS")):
            os.mkdir(os.path.join(self.dir, "TUNNELS"))
        return os.path.join(self.dir, "TUNNELS")

class StructureAssets:

    rcsb_id: str
    paths  : StructureAssetPaths

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id
        self.paths   = StructureAssetPaths(rcsb_id)
    
    def assets_status(self)->dict[AssetClass, bool]:
        _ = {}
        for a in AssetClass:
            _[a]= self.verify_exists(a)
        return _

    def verify_exists(self, asset: AssetClass) -> bool:
        match asset:
            case AssetClass.profile:
                return os.path.exists(self.paths.profile)
            case AssetClass.mmcif:
                return os.path.exists(self.paths.cif)
            case AssetClass.ptc:
                return os.path.exists(self.paths.ptc)
            # case AssetClass.chains:
            #     return os.path.exists(self.paths.chains_dir) and len(os.listdir(self.paths.chains_dir))> 0
            case AssetClass.thumbnail:
                return os.path.exists(self.paths.thumbnail)
            case _:
                raise KeyError("AssetClass {} does not exist.".format(asset))

    async def upsert_cif(self, overwrite: bool = False):
        if not os.path.exists(self.paths.cif):
            await download_unpack_place(self.rcsb_id)
            print("Saved cif file: \t", self.paths.cif)
        else:
            if overwrite:
                await download_unpack_place(self.rcsb_id)
                print("Overwrote cif file: \t", self.paths.cif)

    def upsert_thumbnail(self, overwrite: bool = False) -> bool:
        #TODO: ChimeraX image gen
        ...

    async def upsert_chains(self):
        #TODO: ChimeraX splitchain
        try:
            #TODO: chimerax should be a env variable pointint to the binary, so too should the script
            # Command to run ChimeraX with the script
            command = ['chimerax',   '--nogui', '--offscreen' ,'--cmd' ,'open /home/rtviii/dev/riboxyz/chimerax/chainsplitter.py; open {}; chainsplitter {}; exit'.format(self.rcsb_id, self.rcsb_id)]
            process = subprocess.Popen(
                command,
                stdout = subprocess.PIPE,
                stderr = subprocess.PIPE,
                text   = True
            )
            # Capture output in real-time
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    print("Output:", output.strip())

            # Capture any remaining output and errors
            stdout, stderr = process.communicate()

            if process.returncode == 0:
                print("ChimeraX script executed successfully.")
            else:
                print("ChimeraX script execution failed.")
                print("Error:", stderr)

            return process.returncode, stdout, stderr

        except Exception as e:
            print(f"An error occurred: {str(e)}")
            return None, None, str(e)
        
    async def upsert_ptc(self, overwrite: bool = False):

        etllogger = get_etl_logger()
        asset_ptc_coords_path = os.path.join(self.paths.ptc)

        if os.path.exists(asset_ptc_coords_path) and not overwrite:
            logger.debug(f"PTC coordinates already exist for {self.rcsb_id} and overwrite is set to False" )

        # ress, auth_asym_id = ptc_resdiues_get( self.ro.biopython_structure(), self.ro.profile().rnas )
        # midpoint_coords    = ptc_residues_calculate_midpoint(ress, auth_asym_id)

        # writeout = {
        #     "site_9_residues"      : [(res.get_resname(), res.id[1]) for res in ress],
        #     "LSU_rRNA_auth_asym_id": auth_asym_id,
        #     "midpoint_coordinates" : midpoint_coords,
        #     "nomenclature_table"   : self.ro.nomenclature_table(),
        # }

        # with open(asset_ptc_coords_path, "w") as f:
        #     json.dump(writeout, f)
        #     etllogger.info( f"Saved PTC coordinates for {self.rcsb_id} to {asset_ptc_coords_path}" )
        raise NotImplemented("IMplement and move out of the assests")
        ...

    def write_own_json_profile(self, new_profile: dict, overwrite: bool = False):
        """Update self, basically."""
        if os.path.exists(self.paths.profile) and not overwrite:
            print( "You are about to overwrite {}. Specify `overwrite=True` explicitly.".format( self.paths.profile )
            )
        elif overwrite:
            with open(self.paths.profile, "w") as f:
                json.dump(new_profile, f)
                logger.debug(f"Updated profile for {self.rcsb_id}")

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

    def _verify_dir_exists(self):
        if not os.path.exists(self.paths.dir):
            os.umask(0)
            os.makedirs(self.paths.dir, 0o755)