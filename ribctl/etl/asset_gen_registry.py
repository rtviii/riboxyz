"""New registry system for managing asset generators"""
from pathlib import Path
import functools
from loguru import logger
from typing import Dict, Optional

from .asset_manager import AssetType, AssetGenerator, RibosomeAssetManager

def log_generator(func: AssetGenerator) -> AssetGenerator:
    """Decorator to add consistent logging to generators"""
    @functools.wraps(func)
    async def wrapper(pdb_id: str, output_path: Path, force: bool = False) -> None:
        logger.info(f"Starting {func.__name__} for {pdb_id}")
        try:
            if output_path.exists() and not force:
                logger.info(f"Asset exists at {output_path}, force=False, skipping")
                return
            await func(pdb_id, output_path, force)
            logger.success(f"Generated {func.__name__} for {pdb_id}")
        except Exception as e:
            logger.exception(f"Failed {func.__name__} for {pdb_id}: {str(e)}")
            raise
    return wrapper

class AssetRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager

    def register(self, asset_type: AssetType):
        """Decorator to register a generator function"""
        def decorator(func: AssetGenerator) -> AssetGenerator:
            logged_func = log_generator(func)
            self.manager.register_generator(asset_type, logged_func)
            return logged_func
        return decorator

    async def generate_asset(self, pdb_id: str, asset_type: AssetType, force: bool = False) -> None:
        """Generate a single asset and its dependencies"""
        logger.info(f"Generating {asset_type.name} for {pdb_id}")
        
        asset_def = self.manager.assets[asset_type]
        
        # Handle dependencies first
        if asset_def.dependencies:
            logger.info(f"Processing dependencies for {asset_type.name}")
            for dep in asset_def.dependencies:
                await self.generate_asset(pdb_id, dep, force)

        # Generate the asset
        if not asset_def.generator:
            raise ValueError(f"No generator registered for {asset_type}")

        output_path = self.manager.path_manager.get_asset_path(pdb_id, asset_type)
        await asset_def.generator(pdb_id, output_path, force)

    async def generate_multiple(self, pdb_id: str, asset_types: list[AssetType], force: bool = False) -> None:
        """Generate multiple assets for a structure"""
        for asset_type in asset_types:
            await self.generate_asset(pdb_id, asset_type, force)



# Example Usage:
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