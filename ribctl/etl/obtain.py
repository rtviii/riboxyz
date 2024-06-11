import asyncio
from concurrent import futures
import json
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from typing import Optional
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.etl import AssetFile
from ribctl.etl.ribosome_assets import  Asset, RibosomeAssets
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.schema.types_binding_site import BindingSite
from ribctl.lib.mod_extract_bsites import bsite_ligand, struct_ligand_ids, bsite_extrarbx_polymer, bsite_extrarbx_polymer
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.etl_pipeline import current_rcsb_structs, ReannotationPipeline, rcsb_single_structure_graphql, query_rcsb_api
from concurrent.futures import  Future, ThreadPoolExecutor
from ribctl.logs.loggers import get_etl_logger

def obtain_assets(rcsb_id: str, assetlist:list[Asset], overwrite: bool = False):
    """This should return an array of Futures for acquisition routines for each A  in asset type."""

    rcsb_id = rcsb_id.upper()
    assets  = RibosomeAssets(rcsb_id)
    assets._verify_dir_exists()

    coroutines = []


    if Asset.profile in assetlist:
        coroutines.append(assets.update_profile(overwrite))

    if Asset.cif in assetlist:
        coroutines.append(assets.update_cif(overwrite))

    if Asset.ptc in assetlist:
        coroutines.append(assets.update_ptc(overwrite))

    if Asset.chains in assetlist:
        coroutines.append(...)

    return coroutines

def obtain_assets_threadpool(assetlist, workers: int = 10,  overwrite=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = get_etl_logger()
    unsynced = sorted(current_rcsb_structs())
        
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


# def obtain_all_missing():
#     import asyncio
#     etllogger = get_etl_logger()
#     statuses = AssetFile.status_all()
#     print(statuses)
#     print("got all struct statust", len(statuses))
#     count = 0
#     for i,j in statuses:
#         if j['PTC'] == False:
#             count += 1
#     print("PTCs missing", count)

#     for struct, status in statuses:
#         if not status['PTC'] :
#             try:
#                 asyncio.run(RibosomeAssets(struct)._update_ptc_coordinates())
#             except Exception as e:
#                 etllogger.error(f"Error in {struct} : {e}")