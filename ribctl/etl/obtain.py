import asyncio
from concurrent import futures
import json
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from typing import Optional
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.ribosome_types.types_binding_site import BindingSite
from ribctl.lib.mod_extract_bsites import bsite_ligand, struct_ligand_ids, bsite_extrarbx_polymer, bsite_extrarbx_polymer
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.etl_pipeline import current_rcsb_structs, ReannotationPipeline, rcsb_single_structure_graphql, query_rcsb_api
from concurrent.futures import  Future, ThreadPoolExecutor

import os

from ribctl.logs.loggers import get_etl_logger

async def obtain_assets(rcsb_id: str, assetlist: Assetlist, overwrite: bool = False):
    """Obtain all assets for a given RCSB ID"""

    rcsb_id = rcsb_id.upper()
    assets = RibosomeAssets(rcsb_id)
    assets._verify_dir_exists()

    coroutines = []

    print("Obtaining assets")

    if assetlist.profile:
        print("Obtaining assets:profile")
        coroutines.append(assets._update_json_profile(overwrite))

    if assetlist.cif:
        print("Obtaining assets:cif")
        coroutines.append(assets._update_cif(overwrite))

    if assetlist.ligands:
        print("Obtaining assets:ligands")
        coroutines.append(assets._update_ligands(overwrite))

    if assetlist.ptc_coords:
        print("Obtaining assets:ptc_coords")
        asset_ptc_coords_path = os.path.join(assets._dir_path(),f'{assets.rcsb_id}_PTC_COORDINATES.json')

        if os.path.exists(asset_ptc_coords_path) and not overwrite:
            raise Exception(f'PTC coordinates already exist for {assets.rcsb_id} and overwrite is set to False')

        ress, auth_asym_id = ptc_resdiues_get(assets.biopython_structure(),  assets.profile().rnas)
        midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)

        writeout = {
            "site_9_residues"      : [(res.get_resname(), res.get_segid(), res.full_id) for res in ress],
            "LSU_rRNA_auth_asym_id": auth_asym_id,
            "midpoint_coordinates" : midpoint_coords,
            'nomenclature_table'   : assets._nomenclature_table()
        }

        with open(asset_ptc_coords_path, 'w') as f:
            json.dump(writeout, f)

    if assetlist.cif_modified_and_chains:
        print("Obtaining assets:cif modified and chains")
        coroutines.append(assets.upsert(overwrite))

    await asyncio.gather(*coroutines)

def obtain_assets_threadpool(targets: list[str], assetlist: Assetlist, workers: int = 15, get_all: bool = False, overwrite=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = get_etl_logger()

    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(), targets))
        
    logger.info(f"Found {len(unsynced)} unsynced structures")

    tasks: list[Future] = []
    results = []
    logger.debug("Begun downloading ribosome profiles via RCSB")

    def log_commit_result(rcsb_id: str):
        def _(f: Future):
            if not None == f.exception():
                logger.error(rcsb_id + ":" + f.exception().__str__())
            else:
                logger.debug(rcsb_id + ": processed successfully.")
        return _

    with ThreadPoolExecutor(max_workers=workers) as executor:

        print("Got unsynced structures")
        for rcsb_id in unsynced:
            fut = executor.submit(asyncio.run, obtain_assets(rcsb_id, assetlist, overwrite))
            fut.add_done_callback(log_commit_result(rcsb_id))
            tasks.append(fut)

        for future in futures.as_completed(tasks):
            try:
                results.append(future.result())
            except Exception as e:
                logger.error(future.exception())

    logger.info("Finished syncing with RCSB")
