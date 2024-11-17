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

@registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    profile = await ETLCollector(rcsb_id).process_structure()
    return profile

@registry.register(AssetType.PTC)
async def generate_ptc(rcsb_id: str) -> PTCInfo:
    return PTC_location(rcsb_id)


rcsb_id = '4UG0'
asyncio.run(registry.generate_asset(rcsb_id,AssetType.STRUCTURE_PROFILE))