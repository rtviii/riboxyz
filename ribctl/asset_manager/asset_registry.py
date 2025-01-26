from typing import TypeVar, Callable, Awaitable
from pathlib import Path
import functools
from loguru import logger
from pydantic import BaseModel
from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_raw import RawAssetHandler
from ribctl.etl.etl_collector import ETLCollector
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.schema.types_ribosome import ConstrictionSite, PTCInfo, RibosomeStructure
from .types import AssetType

ModelT = TypeVar("ModelT", bound=BaseModel)

class AssetRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager
        self.path_manager = manager.path_manager
        self.raw_handler = RawAssetHandler(manager.path_manager.base_dir)


    def register(self, asset_type: AssetType):
        def decorator(
            func: Callable[[str], Awaitable[ModelT]]
        ) -> Callable[[str, bool], Awaitable[None]]:
            @functools.wraps(func)
            async def wrapped(rcsb_id: str, overwrite: bool = False) -> None:
                output_path = self.path_manager.get_asset_path(rcsb_id, asset_type)
                try:
                    if output_path.exists() and not overwrite:
                        logger.info(f"Asset exists at {output_path}, skipping")
                        return

                    output_path.parent.mkdir(parents=True, exist_ok=True)

                    # Load only model-based dependencies
                    dependencies = {}
                    if asset_type.dependencies:
                        for dep in asset_type.dependencies:
                            if dep.is_raw_asset:
                                # Skip trying to load raw assets as models
                                continue

                            dep_path = self.path_manager.get_asset_path(rcsb_id, dep)
                            if not dep_path.exists():
                                await self.generate_asset(rcsb_id, dep, overwrite)
                            model_cls = dep.model_type
                            dependencies[dep.value.name] = (
                                model_cls.model_validate_json(dep_path.read_text())
                            )

                    # Generate and save
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
        if asset_type.is_raw_asset:
            await self.raw_handler.handle_asset(rcsb_id, asset_type, force)
        else:
            asset_def = self.manager.assets[asset_type]
            if not asset_def.generator:
                raise ValueError(f"No generator registered for {asset_type}")
            await asset_def.generator(rcsb_id, force)

    async def generate_multiple(
        self, rcsb_id: str, asset_types: list[AssetType], force: bool = False
    ) -> None:
        """Generate multiple assets for a structure"""
        for asset_type in asset_types:
            await self.generate_asset(rcsb_id, asset_type, force)

main_registry = AssetRegistry(RibosomeAssetManager(RIBETL_DATA))

@main_registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    profile = await ETLCollector(rcsb_id).generate_profile(overwrite=False, reclassify=True)
    return profile

@main_registry.register(AssetType.PTC)
async def generate_ptc(rcsb_id: str) -> PTCInfo:
    return PTC_location(rcsb_id)

@main_registry.register(AssetType.CONSTRICTION_SITE)
async def generate_constriction(rcsb_id: str) -> ConstrictionSite:
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
