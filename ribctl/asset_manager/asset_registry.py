from typing import TypeVar, Callable, Awaitable
from pathlib import Path
import functools
from loguru import logger

from pydantic import BaseModel

from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from .asset_types import AssetType

ModelT = TypeVar('ModelT', bound=BaseModel)

class AssetRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager
        self.path_manager = manager.path_manager

    def register(self, asset_type: AssetType) -> Callable[
        [Callable[[str], Awaitable[ModelT]]], 
        Callable[[str, bool], Awaitable[None]]
    ]:
        """Decorator for registering model-based asset generators"""
        def decorator(func: Callable[[str], Awaitable[ModelT]]) -> Callable[[str, bool], Awaitable[None]]:
            @functools.wraps(func)
            async def wrapped(rcsb_id: str, overwrite: bool = False) -> None:
                output_path = self.path_manager.get_asset_path(rcsb_id, asset_type)
                logger.info(f"Starting {func.__name__} for {rcsb_id}")

                try:
                    if output_path.exists() and not overwrite:
                        logger.info(f"Asset exists at {output_path}, skipping")
                        return

                    output_path.parent.mkdir(parents=True, exist_ok=True)

                    # Load dependencies if needed
                    dependencies = {}
                    if asset_type.dependencies:
                        for dep in asset_type.dependencies:
                            dep_path = self.path_manager.get_asset_path(rcsb_id, dep)
                            if not dep_path.exists():
                                await self.generate_asset(rcsb_id, dep, overwrite)
                            model_cls = dep.model_type
                            dependencies[dep.value.name] = model_cls.model_validate_json(dep_path.read_text())

                    # Generate and save the model
                    result = await func(rcsb_id)
                    output_path.write_text(result.model_dump_json())
                    logger.success(f"Generated {asset_type.name} for {rcsb_id}")

                except Exception as e:
                    logger.exception(f"Failed {func.__name__} for {rcsb_id}: {str(e)}")
                    raise

            self.manager.register_generator(asset_type, wrapped)
            return wrapped
        return decorator

    async def generate_asset(self, rcsb_id: str, asset_type: AssetType, force: bool = False) -> None:
        """Generate a single asset and its dependencies"""
        logger.info(f"Generating {asset_type.name} for {rcsb_id}")
        
        asset_def = self.manager.assets[asset_type]
        if not asset_def.generator:
            raise ValueError(f"No generator registered for {asset_type}")

        await asset_def.generator(rcsb_id, force)

    async def generate_multiple(self, rcsb_id: str, asset_types: list[AssetType], force: bool = False) -> None:
        """Generate multiple assets for a structure"""
        for asset_type in asset_types:
            await self.generate_asset(rcsb_id, asset_type, force)