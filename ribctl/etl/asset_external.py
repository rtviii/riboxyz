from pathlib import Path
import asyncio
from loguru import logger
from ribctl.lib.utils import download_unpack_place

class RawAssetHandler:
    """Simple handler for raw file assets like mmcif files"""
    
    def __init__(self, base_dir: Path):
        self.base_dir = base_dir

    async def fetch_mmcif(self, rcsb_id: str, force: bool = False) -> None:
        """Download and save mmCIF file"""
        output_path = self.base_dir / rcsb_id.upper() / f"{rcsb_id}.cif"
        
        if output_path.exists() and not force:
            logger.info(f"MMCIF exists for {rcsb_id}, skipping")
            return
            
        output_path.parent.mkdir(parents=True, exist_ok=True)
        await download_unpack_place(rcsb_id)
        logger.success(f"Downloaded MMCIF for {rcsb_id}")