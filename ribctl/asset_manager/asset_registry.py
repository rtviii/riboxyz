from typing import TypeVar, Callable, Awaitable
from pathlib import Path
import functools
from loguru import logger
from pydantic import BaseModel

from pathlib import Path
import asyncio
from loguru import logger
from typing import Dict, Callable, Awaitable
from functools import partial

from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
# from ribctl.lib.npet.alpah_lib import produce_alpha_contour
# from ribctl.lib.npet.alpha_lib_parallel import produce_alpha_contour
# from ribctl.lib.npet.alpha_lib import produce_alpha_contour
from ribctl.lib.npet.contour_via_poisson_recon import alpha_contour_via_poisson_recon
from ribctl.lib.npet.npet_driver import create_npet_mesh
from ribctl.lib.utils import download_unpack_place
from ribctl.asset_manager.asset_types import AssetType
from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.etl.etl_collector import ETLCollector
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.schema.types_ribosome import (
    ConstrictionSite,
    PTCInfo,
    RibosomeStructure,
)
from .asset_types import AssetType

ModelT = TypeVar("ModelT", bound=BaseModel)


class RawAssetHandler:
    """Handler for raw file assets with extensible asset type matching"""

    def __init__(self):
        self._handlers: Dict[AssetType, Callable[[str, bool], Awaitable[None]]] = {}
        self._register_default_handlers()

    def _register_default_handlers(self) -> None:
        """Register built-in handlers for known raw asset types"""
        self.register_handler(AssetType.MMCIF, self._fetch_mmcif)
        self.register_handler(AssetType.NPET_MESH, npet_mesh_handler)
        self.register_handler(AssetType.ALPHA_SHAPE, alphashape_handler)

    def register_handler(
        self, asset_type: AssetType, handler: Callable[[str, bool], Awaitable[None]]
    ) -> None:
        """Register a new handler for an asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(
                f"Cannot register handler for non-raw asset type: {asset_type}"
            )
        self._handlers[asset_type] = handler

    async def handle_asset(
        self, rcsb_id: str, asset_type: AssetType, force: bool = False
    ) -> None:
        """Generic handler for any registered raw asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(f"Asset type {asset_type} is not a raw asset")

        handler = self._handlers.get(asset_type)
        if not handler:
            raise ValueError(f"No handler registered for raw asset type: {asset_type}")

        await handler(rcsb_id, force)

    async def _fetch_mmcif(self, rcsb_id: str, force: bool = False) -> None:
        """Download and save mmCIF file"""
        output_path = AssetType.MMCIF.get_path(rcsb_id)

        if output_path.exists() and not force:
            logger.info(f"MMCIF exists for {rcsb_id}, skipping")
            return

        output_path.parent.mkdir(parents=True, exist_ok=True)
        await download_unpack_place(rcsb_id)
        logger.success(f"Downloaded MMCIF for {rcsb_id}")

    # Example of how to add another handler:
    # async def _fetch_npet_mesh(self, rcsb_id: str, force: bool = False) -> None:
    #     """Download and save NPET mesh file"""
    #     output_path = self.base_dir / rcsb_id.upper() / "TUNNELS" / f"{rcsb_id}_NPET_MESH.ply"
    #
    #     if output_path.exists() and not force:
    #         logger.info(f"NPET mesh exists for {rcsb_id}, skipping")
    #         return
    #
    #     output_path.parent.mkdir(parents=True, exist_ok=True)
    #     # Add actual download/generation logic here
    #     logger.success(f"Generated NPET mesh for {rcsb_id}")


class AssetRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager
        self.raw_handler = RawAssetHandler()

    def register(self, asset_type: AssetType):
        def decorator(
            func: Callable[[str], Awaitable[ModelT]]
        ) -> Callable[[str, bool], Awaitable[None]]:
            @functools.wraps(func)
            async def wrapped(rcsb_id: str, overwrite: bool = False) -> None:
                output_path = asset_type.get_path(rcsb_id)
                try:
                    if output_path.exists() and not overwrite:
                        logger.info(f"Asset exists at {output_path}, skipping")
                        return

                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    result = await func(rcsb_id)
                    output_path.write_text(result.model_dump_json())
                    logger.success(f"Generated {asset_type.name} for {rcsb_id}")

                except Exception as e:
                    logger.exception(f"Failed {func.__name__} for {rcsb_id}: {str(e)}")
                    raise

            self.manager.register_generator(asset_type, wrapped)
            return wrapped

        return decorator

    async def generate_asset(
        self, rcsb_id: str, asset_type: AssetType, force: bool = False
    ) -> None:
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


async def npet_mesh_handler(rcsb_id: str, force: bool) -> None:
    create_npet_mesh(rcsb_id)

async def alphashape_handler(rcsb_id: str, force: bool) -> None:
    alpha_contour_via_poisson_recon(rcsb_id)

main_registry = AssetRegistry(RibosomeAssetManager(RIBETL_DATA))

@main_registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    profile = await ETLCollector(rcsb_id).generate_profile(
        overwrite=False, reclassify=True
    )
    return profile


@main_registry.register(AssetType.PTC)
async def generate_ptc(rcsb_id: str) -> PTCInfo:
    return PTC_location(rcsb_id)


@main_registry.register(AssetType.CONSTRICTION_SITE)
async def generate_constriction(rcsb_id: str) -> ConstrictionSite:
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
