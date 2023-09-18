import asyncio
import json
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from typing import Optional
from api.logs.loggers import get_updates_logger
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.lib.tunnel.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.types.types_binding_site import BindingSite
from ribctl.lib.types.types_poly_nonpoly_ligand import PolymericFactorClass, RNAClass
from ribctl.lib.mod_extract_bsites import bsite_nonpolymeric_ligand, struct_ligand_ids, struct_polymeric_factor_ids, bsite_polymeric_factor, bsite_polymeric_factor
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.etl_pipeline import current_rcsb_structs, ReannotationPipeline, rcsb_single_structure_graphql, query_rcsb_api
from ribctl.lib.mod_render_thumbnail import render_thumbnail
from ribctl.lib.utils import download_unpack_place, open_structure
from ribctl.lib.types.types_ribosome import RNA, PolymerClass, PolymericFactor, Protein, ProteinClass, RibosomeStructure
from ribctl import RIBETL_DATA
from pydantic import BaseModel, parse_obj_as
from concurrent.futures import ALL_COMPLETED, Future, ProcessPoolExecutor, ThreadPoolExecutor, wait
import os

class Assetlist(BaseModel)   : 
      profile                : Optional[bool]
      ptc_coords             : Optional[bool]
      cif                    : Optional[bool]
      cif_modified_and_chains: Optional[bool]
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
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    def _cif_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}.cif"

    def _cif_modified_filepath(self):
        self._envcheck()
        return f"{self._dir_path()}/{self.rcsb_id}_modified.cif"

    def _ptc_residues(self) -> dict[str, dict[str, list[float]]]:
        PTC_RESIDUES_PATH = os.path.join(RIBETL_DATA, self.rcsb_id, "{}_PTC_COORDINATES.json".format(self.rcsb_id))
        with open(PTC_RESIDUES_PATH, 'r') as infile:
            return json.load(infile)

    def ____nomenclature_v2(self) -> dict[str, ProteinClass]:
        with open("/home/rxz/dev/docker_ribxz/api/ribctl/assets/nomenclaturev2/{}.json".format(self.rcsb_id.upper()), 'r') as infile:
            return json.load(infile)

    def _json_profile_filepath(self):
        self._envcheck()
        return os.path.join(self._dir_path(),f"{self.rcsb_id}.json")

    def profile(self) -> RibosomeStructure:
        with open(self._json_profile_filepath(), "r") as f:
            return RibosomeStructure.parse_obj(json.load(f))

    def nomenclature_table(self) -> dict[str, ProteinClass]:
        prof = self.profile()
        m    = {}

        for prot in prof.proteins:
            m[prot.auth_asym_id] = prot.nomenclature

        if prof.rnas!=None:
            for rna in prof.rnas:
                m[rna.auth_asym_id] = rna.nomenclature
        return m
                

        


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

    @staticmethod
    def list_all_structs():
        return os.listdir(RIBETL_DATA)

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Getters =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def get_taxids(self) -> tuple[list[int], list[int]]:
        p = self.profile()
        return (p.src_organism_ids, p.host_organism_ids)

    def get_struct_and_profile(self) -> tuple[Structure, RibosomeStructure]:
        return self.biopython_structure(), self.profile()

    def get_chain_by_polymer_class(self, poly_class: PolymerClass | PolymericFactorClass, assembly: int = 0) -> PolymericFactor | RNA | Protein | None:
        profile = self.profile()
        for prot in profile.proteins:
            if poly_class in prot.nomenclature and prot.assembly_id == assembly:
                return prot
        if profile.rnas is not None:
            for rna in profile.rnas:
                if poly_class in rna.nomenclature and rna.assembly_id == assembly:
                    return rna
        if profile.polymeric_factors is not None:
            for polyf in profile.polymeric_factors:
                if poly_class in polyf.nomenclature and polyf.assembly_id == assembly:
                    return polyf
        return None

    def get_chain_by_auth_asym_id(self, auth_asym_id: str) -> tuple[
            RNA | Protein | PolymericFactor | None,
            typing.Literal["RNA", "Protein", "PolymericFactor"] | None]:

        profile = self.profile()
        for chain in profile.proteins:
            if chain.auth_asym_id == auth_asym_id:
                return (chain, "Protein")

        if profile.rnas is not None:
            for chain in profile.rnas:
                if chain.auth_asym_id == auth_asym_id:
                    return (chain, "RNA")

        if profile.polymeric_factors is not None:
            for chain in profile.polymeric_factors:
                if chain.auth_asym_id == auth_asym_id:
                    return (chain, "PolymericFactor")

        return (None, None)

    def get_rna_by_nomclass(self, class_: RNAClass, assembly: int = 0) -> RNA | None:
        """@assembly here stands to specify which of the two or more models the rna comes from
        in the case that a structure contains multiple models (ex. 4V4Q XRAY)"""

        profile = self.profile()

        if profile.rnas == None:
            return None

        for rna in profile.rnas:
            if class_ in rna.nomenclature and rna.assembly_id == assembly:
                return rna

    def get_prot_by_nomclass(self, class_: ProteinClass, assembly: int = 0) -> Protein | None:
        _auth_asym_id = {
            v: k for k, v in self.____nomenclature_v2().items()}.get(class_, None)
        if _auth_asym_id != None:
            return self.get_chain_by_auth_asym_id(_auth_asym_id)[0]

        for prot in self.profile().proteins:
            if class_ in prot.nomenclature and prot.assembly_id == assembly:
                return prot
        else:
            return None

    def get_LSU_rRNA(self, assembly: int = 0) -> RNA:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """

        rna = self.get_rna_by_nomclass("23SrRNA", assembly)
        if rna == None:
            rna = self.get_rna_by_nomclass("25SrRNA", assembly)
        if rna == None:
            rna = self.get_rna_by_nomclass("28SrRNA", assembly)
        if rna == None:
            raise Exception("No LSU rRNA found in structure")
        else:
            return rna

    def biopython_get_chain(self, auth_asym_id: str) -> Chain:
        return self.biopython_structure().child_dict[0].child_dict[auth_asym_id]

    @staticmethod
    def biopython_chain_get_seq(struct: Structure, auth_asym_id: str, protein_rna: typing.Literal["protein", "rna"], sanitized: bool = False) -> str:

        chain3d = struct.child_dict[0].child_dict[auth_asym_id]
        ress = chain3d.child_list

        # if sanitized == True:
        #     print("sanitized")
        #     ress = [*filter(lambda r: r.get_resname()
        #                     in ["A", "C", "G", "U", "PSU"], ress)]

        #     for _r in ress:
        #         if _r.get_resname() == "PSU":
        #             _r.resname = "U"
        seq = ''
        for i in ress:
            if i.resname not in AMINO_ACIDS_3_TO_1_CODE.keys():
                print("Unknown residue:", i.resname)
                continue

            if protein_rna == 'rna':
                seq += i.resname
            else:
                seq += AMINO_ACIDS_3_TO_1_CODE[i.resname]

        # return reduce(lambda x, y: x + y.resname if protein_rna == 'rna' else AMINO_ACIDS_3_TO_1_CODE[y.resname], ress, '')
        return seq

    # ※ -=-=-=-=-=-=-=-=-=-=-=-=-=-=-= Verification =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

    def _verify_dir_exists(self):
        if not os.path.exists(self._dir_path()):
            os.umask(0)
            os.makedirs(self._dir_path(), 0o777)

    async def _verify_cif(self, overwrite: bool = False) -> bool:
        if not os.path.exists(self._cif_filepath()):
            await download_unpack_place(self.rcsb_id)
            print("Saved structure file:\t", self._cif_filepath())
            return True
        else:
            if overwrite:
                await download_unpack_place(self.rcsb_id)
                print("Saved structure file:\t", self._cif_filepath())
                return True
            else:
                return False

    async def _verify_cif_modified_and_chains(self, overwrite: bool = False) -> bool:
        if not os.path.isdir(self.chains_dir()):
            os.makedirs(self.chains_dir())
            await split_rename(self.rcsb_id)
            return True
        else:
            if overwrite:
                await split_rename(self.rcsb_id)
                return True
            else:
                return False

    async def _verify_json_profile(self, overwrite: bool = False) -> bool:
        self._verify_dir_exists()

        if not os.path.exists(self._json_profile_filepath()):
            ribosome = ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(self.rcsb_id.upper()))).process_structure()
            if not parse_obj_as(RibosomeStructure, ribosome):
                raise Exception("Invalid ribosome structure profile.")

            self.save_json_profile(
                self._json_profile_filepath(), ribosome.dict())
            print(f"Wrote structure profile:\t{self._json_profile_filepath()}")
            return True;
        else:
            if overwrite:
                ribosome = ReannotationPipeline(query_rcsb_api(rcsb_single_structure_graphql(self.rcsb_id.upper()))).process_structure()
                if not parse_obj_as(RibosomeStructure, ribosome):
                    raise Exception("Invalid ribosome structure profile.")
                self.save_json_profile(
                    self._json_profile_filepath(), ribosome.dict())
                print(
                    f"Wrote structure profile:\t{self._json_profile_filepath()}")
                return True
            else:
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

        ligands = struct_ligand_ids(self.rcsb_id, self.profile())
        polymeric_factors = struct_polymeric_factor_ids(self.profile())
        all_verified_flag = True

        for ligand_chemid in ligands:
            if not os.path.exists(BindingSite.path_nonpoly_ligand(self.rcsb_id, ligand_chemid)):
                all_verified_flag = False
                bsite = bsite_nonpolymeric_ligand(
                    ligand_chemid, self.biopython_structure())
                bsite.save(bsite.path_nonpoly_ligand(
                    self.rcsb_id, ligand_chemid))
            else:
                if overwrite:
                    bsite = bsite_nonpolymeric_ligand(
                        ligand_chemid, self.biopython_structure())
                    bsite.save(bsite.path_nonpoly_ligand(
                        self.rcsb_id, ligand_chemid))
                else:
                    ...

        if polymeric_factors is not None:
            for poly in polymeric_factors:
                if not os.path.exists(BindingSite.path_poly_factor(self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id)):
                    all_verified_flag = False
                    bsite = bsite_polymeric_factor(
                        poly.auth_asym_id, self.biopython_structure())
                    bsite.save(bsite.path_poly_factor(
                        self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id))
                else:
                    if overwrite:
                        bsite = bsite_polymeric_factor(
                            poly.auth_asym_id, self.biopython_structure())
                        bsite.save(bsite.path_poly_factor(
                            self.rcsb_id, poly.nomenclature[0], poly.auth_asym_id))
                    else:
                        ...

        return all_verified_flag

# ※ Mass process methods.

async def obtain_assets(rcsb_id: str, assetlist: Assetlist, overwrite: bool = False):
    """Obtain all assets for a given RCSB ID"""

    assets = RibosomeAssets(rcsb_id)
    assets._verify_dir_exists()

    coroutines = []

    if assetlist.profile:
        coroutines.append(assets._verify_json_profile(overwrite))

    if assetlist.cif:
        coroutines.append(assets._verify_cif(overwrite))

    if assetlist.factors_and_ligands:
        coroutines.append(assets._verify_ligads_and_ligandlike_polys(overwrite))

    if assetlist.ptc_coords:
            ress, auth_asym_id = ptc_resdiues_get(RCSB_ID)
            midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)

            residue_labels = [(res.get_resname(), res.id[1]) for res in ress]
            print(residue_labels)

            writeout = {
                "site_9_residues"      : [(res.get_resname(), res.get_segid()) for res in ress],
                "LSU_rRNA_auth_asym_id": auth_asym_id,
                "midpoint_coordinates" : midpoint_coords,
                'nomenclature_table'   : assets.nomenclature_table()
            }

            with open(os.path.join(assets._dir_path(),f'{assets.rcsb_id}_PTC_COORDINATES.json'), 'w') as f:
                json.dump(writeout, f)


    if assetlist.cif_modified_and_chains:
        coroutines.append(assets._verify_cif_modified_and_chains(overwrite))

    await asyncio.gather(*coroutines)

def obtain_assets_threadpool(targets: list[str], assetlist: Assetlist, workers: int = 5, get_all: bool = False, overwrite=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = get_updates_logger()
    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(), targets))

    futures: list[Future] = []
    logger.info("Begun downloading ribosome profiles via RCSB")

    def log_commit_result(rcsb_id: str):
        def _(f: Future):
            if not None == f.exception():
                logger.error(rcsb_id + ":" + f.exception().__str__())
            else:
                logger.info(rcsb_id + ": processed successfully.")
        return _

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for rcsb_id in unsynced:
            fut = executor.submit(asyncio.run, obtain_assets(rcsb_id, assetlist, overwrite))
            fut.add_done_callback(log_commit_result(rcsb_id))
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)
    logger.info("Finished syncing with RCSB")

def obtain_assets_processpool(targets: list[str], assetlist: Assetlist, workers: int = 5, get_all: bool = False, overwrite=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = get_updates_logger()
    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(), targets))

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
            fut = executor.submit(asyncio.run, obtain_assets(
                rcsb_id, assetlist, overwrite))
            fut.add_done_callback(log_commit_result(rcsb_id))
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)
    logger.info("Finished syncing with RCSB")
