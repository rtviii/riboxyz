import asyncio
from typing import Coroutine
from ribctl.etl.etl_assets_ops import AssetClass, Assets, RibosomeOps, Structure
from ribctl.etl.etl_pipeline import ( ReannotationPipeline )
import asyncio


# This should be in the RibosomeAssets module
def asset_routines(
    rcsb_id: str, assetlist: list[AssetClass], overwrite: bool = False
) -> list[Coroutine]:
    """This should return an array of Futures for acquisition routines for each A  in asset type."""

    rcsb_id = rcsb_id.upper()
    RO = RibosomeOps(rcsb_id)
    RO._verify_dir_exists()
    RA = Assets(rcsb_id)

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

