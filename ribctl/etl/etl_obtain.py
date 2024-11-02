import asyncio
from typing import Coroutine
from ribctl.etl.assets_structure import AssetClass
from ribctl.ribosome_ops import  StructureAssets, RibosomeOps, Structure
from ribctl.etl.etl_collector import ( ETLCollector )
import asyncio


# This should be in the RibosomeAssets module
def asset_routines(
    rcsb_id: str, assetlist: list[AssetClass], overwrite: bool = False,
    reclassify = False
) -> list[Coroutine]:
    """This should return an array of Futures for acquisition routines for each Asset  in asset type."""

    rcsb_id = rcsb_id.upper()
    RO      = RibosomeOps(rcsb_id)
    RO.assets._verify_dir_exists()
    RA = StructureAssets(rcsb_id)

    coroutines = []
    if AssetClass.profile in assetlist:
        coroutines.append( ETLCollector(rcsb_id).process_structure(overwrite, reclassify) )

    if AssetClass.mmcif in assetlist:
        coroutines.append(RA.upsert_cif(overwrite))

    if AssetClass.ptc in assetlist:
        coroutines.append(RA.upsert_ptc(overwrite))

    return coroutines

async def execute_asset_task_pool(tasks):
    asyncio.gather(*tasks)

