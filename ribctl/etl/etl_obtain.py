import asyncio
from concurrent import futures
from typing import Coroutine, Optional
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.etl.etl_ribosome_ops import AssetClass, Structure
from ribctl.etl.etl_pipeline import (
    current_rcsb_structs,
    ReannotationPipeline,
    rcsb_single_structure_graphql,
    query_rcsb_api,
)
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor
from ribctl.logs.loggers import get_etl_logger
import asyncio


# This should be in the RibosomeAssets module
def asset_routines(
    rcsb_id: str, assetlist: list[AssetClass], overwrite: bool = False
) -> list[Coroutine]:
    """This should return an array of Futures for acquisition routines for each A  in asset type."""

    rcsb_id = rcsb_id.upper()
    RA = Structure(rcsb_id)
    RA._verify_dir_exists()

    coroutines = []

    if AssetClass.profile in assetlist:
        coroutines.append(
            ReannotationPipeline( ReannotationPipeline.rcsb_request_struct(rcsb_id) ).process_structure(overwrite) )

    if AssetClass.cif in assetlist:
        coroutines.append(RA.upsert_cif(overwrite))

    if AssetClass.ptc in assetlist:
        coroutines.append(RA.upsert_ptc(overwrite))

    if AssetClass.chains in assetlist:
        coroutines.append(...)  # todo: chimerax split chains (get 1.8 build)

    return coroutines


async def execute_asset_task_pool(tasks):
    asyncio.gather(*tasks)


def obtain_asssets_threadpool(assetlist, workers: int = 10, overwrite=False):
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
            fut = executor.submit(
                asyncio.run, asset_routines(rcsb_id, assetlist, overwrite)
            )
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
