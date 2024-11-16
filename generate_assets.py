
"""
# Set up manager and registry
manager = RibosomeAssetManager(RIBETL_DATA)
registry = AssetRegistry(manager)

# Register generators with decorator
@registry.register(AssetType.MMCIF)
async def generate_mmcif(pdb_id: str, output_path: Path, force: bool = False) -> None:
    # Your mmCIF generation logic
    ...

# Generate assets
await registry.generate_multiple("1J5E", [AssetType.MMCIF], force=False)
"""

import asyncio
from pathlib import Path
from ribctl import RIBETL_DATA
from ribctl.etl.asset_gen_registry import AssetRegistry
from ribctl.etl.asset_manager import AssetType, RibosomeAssetManager
from ribctl.etl.etl_collector import ETLCollector


manager  = RibosomeAssetManager(RIBETL_DATA)
registry = AssetRegistry(manager)


rcsb_id = "3J7Z"


@registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(pdb_id: str, output_path: Path, force: bool = False) -> None:
    profile = await ETLCollector(rcsb_id=rcsb_id).process_structure(overwrite=force)
    print(profile)
    


asyncio.run( registry.generate_asset(rcsb_id, AssetType.STRUCTURE_PROFILE, force=True))