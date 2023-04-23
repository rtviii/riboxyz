import asyncio
import json
from typing import Optional, Tuple
from api.ribctl.lib.types.types_binding_site import BindingSite
from api.ribctl.lib.types.types_poly_nonpoly_ligand import RNAClass
from logs.loggers import updates_logger
from ribctl.lib.mod_extract_bsites import bsite_nonpolymeric_ligand, struct_ligand_ids, struct_polymeric_factor_ids, bsite_polymeric_factor, bsite_polymeric_factor
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.struct_rcsb_api import current_rcsb_structs, process_pdb_record
from ribctl.lib.mod_render_thumbnail import render_thumbnail
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.types.types_ribosome import ProteinClass, RibosomeStructure
from pydantic import BaseModel, parse_obj_as
import os
from concurrent.futures import ALL_COMPLETED, Future, ProcessPoolExecutor, ThreadPoolExecutor, wait


RIBETL_DATA = str(os.environ.get("RIBETL_DATA"))


class Assetlist(BaseModel): 
      profile                : Optional[bool]
      structure              : Optional[bool]
      structure_modified     : Optional[bool]
      chains_and_modified_cif: Optional[bool]
      factors_and_ligands    : Optional[bool]
      png_thumbnail          : Optional[bool]


class RibosomeAssets():
    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()

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

    def profile(self) -> RibosomeStructure:
        with open(self._json_profile_filepath(), "r") as f:
            return RibosomeStructure.parse_obj(json.load(f))

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

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Getters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def get_seq_rna_by_nomclass(self, class_: RNAClass, canonical:bool=True) -> Tuple[ str, str] | Tuple[None,None]:
        profile = self.profile()

        if profile.rnas == None:
            return ( None,None )

        for rna in profile.rnas:
            if class_ in rna.nomenclature :
                if canonical:
                    return ( rna.entity_poly_seq_one_letter_code_can, rna.auth_asym_id )
                else:
                    return ( rna.entity_poly_seq_one_letter_code, rna.auth_asym_id )
        return ( None,None )

    def get_seq_prot_by_nomclass(self, class_: ProteinClass, canonical:bool=True)-> Tuple[ str,str] | Tuple[None,None] :
        for prot in self.profile().proteins:
            if class_ in prot.nomenclature :
                if canonical:
                    return ( prot.entity_poly_seq_one_letter_code_can, prot.auth_asym_id )
                else:
                    return ( prot.entity_poly_seq_one_letter_code, prot.auth_asym_id )
        else:
            return ( None,None )

    def get_LSU_rRNA(self, get_canonical:bool=True)->Tuple[str,str,str]:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """
        annotated_cifpath = self._cif_modified_filepath()
        rna_type          = ""

        [seq, auth_asym_id] = self.get_seq_rna_by_nomclass("23SrRNA", get_canonical)
        rna_type = "23SrRNA"

        if seq == None or auth_asym_id == None:
            [chain_id, strand_target] = self.get_seq_rna_by_nomclass(
                "25SrRNA",
                get_canonical
            )

            rna_type = "25SrRNA"

        if seq == None or auth_asym_id == None:
            [chain_id, strand_target] = self.get_seq_rna_by_nomclass(
                "28SrRNA",
                get_canonical)
            rna_type = "28SrRNA"

        if seq == None or auth_asym_id == None:
            raise Exception("No LSU rRNA found in structure")

        return (seq,auth_asym_id, rna_type)



    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Verification =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def _verify_dir_exists(self):
        if not os.path.exists(self._dir_path()):
            os.umask(0)
            os.makedirs(self._dir_path(), 0o777)

    async def _verify_cif(self, overwrite: bool = False) -> bool:
        if overwrite:
            download_unpack_place(self.rcsb_id)
            print("Saved structure file:\t", self._cif_filepath())
            return True
        else:
            if os.path.exists(self._cif_filepath()):
                return True
            else:
                return False

    async def _verify_cif_modified_and_chains(self, overwrite: bool = False) -> bool:
        if overwrite:
            split_rename(self.rcsb_id)
            return True
        else:
            return os.path.exists(self._cif_modified_filepath()) and os.path.isdir(self.chains_dir())

    async def _verify_json_profile(self, overwrite: bool = False) -> bool:
        if overwrite:
            ribosome = process_pdb_record(self.rcsb_id)

            if not parse_obj_as(RibosomeStructure, ribosome):
                raise Exception("Invalid ribosome structure profile.")

            self._verify_dir_exists()
            self.save_json_profile(self._json_profile_filepath(), ribosome.dict())
            print(f"Saved structure profile:\t{self._json_profile_filepath()}")
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

    async def _verify_chains_dir(self):
        split_rename(self.rcsb_id)

    async def _verify_ligads_and_ligandlike_polys(self, overwrite: bool = False):

        # def ligand_path(chem_id): return os.path.join(self._dir_path(), f"polymer_{chem_id.upper()}.json")
        # def poly_factor_path(auth_asym_id): return os.path.join(self._dir_path(), f"polymer_{auth_asym_id.upper()}.json")

        ligands           = struct_ligand_ids(self.rcsb_id, self.profile())
        polymeric_factors = struct_polymeric_factor_ids(self.profile())
        all_verified_flag = True

        for ligand_chemid in ligands:
            if not os.path.exists(BindingSite.path_nonpoly_ligand(self.rcsb_id,ligand_chemid)):
                all_verified_flag = False
                bsite = bsite_nonpolymeric_ligand(ligand_chemid, self.biopython_structure())
                bsite.save(bsite.path_nonpoly_ligand(self.rcsb_id, ligand_chemid))
            else:
                if overwrite:
                    bsite = bsite_nonpolymeric_ligand(ligand_chemid, self.biopython_structure())
                    bsite.save(bsite.path_nonpoly_ligand(self.rcsb_id, ligand_chemid))
                else:
                    ...


                    
        if polymeric_factors is not None:
            for poly in polymeric_factors:
                if not os.path.exists(BindingSite.path_poly_factor(self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id)):
                    all_verified_flag = False
                    bsite = bsite_polymeric_factor(poly.auth_asym_id, self.biopython_structure())
                    bsite.save(bsite.path_poly_factor(self.rcsb_id, poly.nomenclature[0],poly.auth_asym_id))
                else:
                    if overwrite:
                        bsite = bsite_polymeric_factor(poly.auth_asym_id, self.biopython_structure())
                        bsite.save(bsite.path_poly_factor(self.rcsb_id, poly.nomenclature[0],poly.auth_asym_id))
                    else:
                        ...


        

        return all_verified_flag



# ※ Mass process methods.


async def obtain_assets(rcsb_id: str ,assetlist:Assetlist ,overwrite: bool = False):
    """Obtain all assets for a given RCSB ID"""

    assets = RibosomeAssets(rcsb_id)
    assets._verify_dir_exists()
    
    coroutines = []

    if assetlist.profile:
        print("p")
        coroutines.append(assets._verify_json_profile(overwrite))

    if assetlist.structure:
        print("s")
        coroutines.append(assets._verify_cif(overwrite))

    if assetlist.factors_and_ligands:
        print("fl")
        coroutines.append(assets._verify_ligads_and_ligandlike_polys(overwrite))

    if assetlist.chains_and_modified_cif:
        print("modified")
        coroutines.append(assets._verify_cif_modified_and_chains(overwrite))

    await asyncio.gather(*coroutines)

def obtain_assets_threadpool(targets:list[str], assetlist:Assetlist, workers: int=5, get_all:bool=False, replace=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = updates_logger
    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(),targets))

    futures: list[Future] = []
    logger.info("Begun downloading ribosome profiles via RCSB")

    def log_commit_result(rcsb_id: str):
        def _(f: Future):
            if not None == f.exception():
                logger.error(rcsb_id + ":" + f.exception().__str__())
            else:
                logger.info(rcsb_id + ": processed profile successfully.")
        return _

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for rcsb_id in unsynced:
            fut = executor.submit(asyncio.run, obtain_assets(rcsb_id,assetlist,replace))
            fut.add_done_callback(log_commit_result(rcsb_id))
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)
    logger.info("Finished syncing with RCSB")

def obtain_assets_processpool(targets:list[str], assetlist:Assetlist, workers: int=5, get_all:bool=False, replace=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = updates_logger
    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(),targets))

    futures: list[Future] = []
    logger.info("Begun downloading ribosome profiles via RCSB")

    def log_commit_result(rcsb_id: str):
        def _(f: Future):
            if not None == f.exception():
                logger.error(rcsb_id + ":" + f.exception().__str__())
            else:
                logger.info(rcsb_id + ": processed profile successfully.")
        return _

    with ProcessPoolExecutor(max_workers=workers) as executor:
        for rcsb_id in unsynced:
            fut = executor.submit(asyncio.run, obtain_assets(rcsb_id,assetlist,replace))
            fut.add_done_callback(log_commit_result(rcsb_id))
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)
    logger.info("Finished syncing with RCSB")
