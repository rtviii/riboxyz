from pathlib import Path
import asyncio
from loguru import logger
from typing import Dict, Callable, Awaitable
from functools import partial

from ribctl.lib.utils import download_unpack_place
from ribctl.asset_manager.asset_types import AssetType

class RawAssetHandler:
    """Handler for raw file assets with extensible asset type matching"""
    
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir
        self._handlers: Dict[AssetType, Callable[[str, bool], Awaitable[None]]] = {}
        self._register_default_handlers()

    def _register_default_handlers(self) -> None:
        """Register built-in handlers for known raw asset types"""
        self.register_handler(AssetType.MMCIF, self._fetch_mmcif)
        # Add other default handlers:
        # self.register_handler(AssetType.NPET_MESH, self._fetch_npet_mesh)
        # self.register_handler(AssetType.THUMBNAIL, self._fetch_thumbnail)

    def register_handler(
        self, 
        asset_type: AssetType, 
        handler: Callable[[str, bool], Awaitable[None]]
    ) -> None:
        """Register a new handler for an asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(f"Cannot register handler for non-raw asset type: {asset_type}")
        self._handlers[asset_type] = handler

    async def handle_asset(self, rcsb_id: str, asset_type: AssetType, force: bool = False) -> None:
        """Generic handler for any registered raw asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(f"Asset type {asset_type} is not a raw asset")
            
        handler = self._handlers.get(asset_type)
        if not handler:
            raise ValueError(f"No handler registered for raw asset type: {asset_type}")
            
        await handler(rcsb_id, force)

    async def _fetch_mmcif(self, rcsb_id: str, force: bool = False) -> None:
        """Download and save mmCIF file"""
        output_path = self.base_dir / rcsb_id.upper() / f"{rcsb_id}.cif"
        
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