import asyncio
from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_external import RawAssetHandler
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_registry import AssetRegistry
from ribctl.asset_manager.asset_types import AssetType
from ribctl.etl.etl_collector import ETLCollector
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.schema.types_ribosome import PTCInfo, RibosomeStructure
from ribctl.lib.utils import download_unpack_place


manager  = RibosomeAssetManager(RIBETL_DATA)
registry = AssetRegistry(manager)
RCSB_ID  = '3J7Z'

@registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    profile = await ETLCollector(rcsb_id).process_structure()
    return profile

@registry.register(AssetType.PTC)
async def generate_ptc(rcsb_id: str) -> PTCInfo:
    return PTC_location(rcsb_id)

# raw_handler = RawAssetHandler(RIBETL_DATA)
# asyncio.run(  raw_handler.fetch_mmcif(RCSB_ID, force=True) )
# asyncio.run(registry.generate_asset(RCSB_ID, AssetType.STRUCTURE_PROFILE , True))
asyncio.run(registry.generate_asset(RCSB_ID, AssetType.PTC ))