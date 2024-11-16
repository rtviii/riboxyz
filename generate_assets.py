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
import os
from loguru import logger
from ribctl import RIBETL_DATA
from ribctl.etl.asset_manager import (
    AssetPathManager,
    AssetType,
    RibosomeAssetManager,
    AssetGenRegistry,
)
from ribctl.etl.etl_collector import ETLCollector
from ribctl.ribosome_ops import RibosomeOps

registry = AssetGenRegistry(RibosomeAssetManager(RIBETL_DATA))

@registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str, overwrite: bool = False) -> None:
    RA = RibosomeOps(rcsb_id)
    if os.path.isfile(RA.assets.paths.profile):
        logger.debug("Profile already exists for {}.".format(rcsb_id))
        if not overwrite:
            return
    profile = await ETLCollector(rcsb_id).process_structure(overwrite=overwrite)
    RA.assets.write_own_json_profile(profile.model_dump(), overwrite=overwrite)



asyncio.run(registry.generate_asset("3J7Z", AssetType.STRUCTURE_PROFILE, force=True))
