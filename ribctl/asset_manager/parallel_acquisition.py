# parallel_acquisition.py
import asyncio
from dataclasses import dataclass
from typing import List, Dict, Set, Optional
from pathlib import Path
from loguru import logger
from concurrent.futures import ProcessPoolExecutor

from ribctl.asset_manager.asset_raw import RawAssetHandler
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_registry import AssetRegistry
from ribctl.asset_manager.asset_types import AssetType
from ribctl.etl.etl_collector import ETLCollector
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.schema.types_ribosome import PTCInfo, RibosomeStructure

@dataclass
class AcquisitionResult:
    rcsb_id: str
    asset_type_name: str  # Store enum name instead of enum object
    success: bool
    error: Optional[str] = None

def process_chunk(

    base_dir                 : str,          # Pass as string instead of Path
    rcsb_ids                 : List[str],
    asset_type_names         : List[str],    # Pass enum names instead of enum objects
    force                    : bool = False,
    max_concurrent_structures: int = 4,
    max_concurrent_assets    : int = 3

) -> Dict[str, List[AcquisitionResult]]:
    """Process a chunk of structures in a separate process"""
    
    async def _process_chunk_async() -> Dict[str, List[AcquisitionResult]]:
        # Reconstruct objects inside the process
        base_path = Path(base_dir)
        asset_types = [AssetType[name] for name in asset_type_names]
        
        raw_handler = RawAssetHandler(base_path)
        manager     = RibosomeAssetManager(base_path)
        registry    = AssetRegistry(manager)
        
        # Register your asset generators here
        @registry.register(AssetType.STRUCTURE_PROFILE)
        async def generate_profile(rcsb_id: str) -> 'RibosomeStructure':
            profile = await ETLCollector(rcsb_id).process_structure()
            return profile

        @registry.register(AssetType.PTC)
        async def generate_ptc(rcsb_id: str) -> 'PTCInfo':
            return PTC_location(rcsb_id)

        async def acquire_asset(rcsb_id: str, asset_type: AssetType) -> AcquisitionResult:
            try:
                if asset_type.is_raw_asset:
                    await raw_handler.handle_asset(rcsb_id, asset_type, force)
                else:
                    await registry.generate_asset(rcsb_id, asset_type, force)
                return AcquisitionResult(rcsb_id, asset_type.name, True)
            except Exception as e:
                logger.exception(f"Failed to acquire {asset_type.name} for {rcsb_id}")
                return AcquisitionResult(rcsb_id, asset_type.name, False, str(e))


        # Process structures with concurrency limits
        structure_sem = asyncio.Semaphore(max_concurrent_structures)
        asset_sem = asyncio.Semaphore(max_concurrent_assets)
        
        async def process_structure(rcsb_id: str) -> tuple[str, List[AcquisitionResult]]:
            async with structure_sem:
                # Group assets by dependency level
                dependency_levels: Dict[int, Set[AssetType]] = {0: set()}
                for asset_type in asset_types:
                    level = len(asset_type.dependencies) if asset_type.dependencies else 0
                    if level not in dependency_levels:
                        dependency_levels[level] = set()
                    dependency_levels[level].add(asset_type)

                all_results = []
                # Process each dependency level
                for level in sorted(dependency_levels.keys()):
                    current_assets = dependency_levels[level]
                    
                    async def acquire_with_semaphore(asset_type: AssetType):
                        async with asset_sem:
                            return await acquire_asset(rcsb_id, asset_type)
                    
                    tasks = [acquire_with_semaphore(asset_type) 
                            for asset_type in current_assets]
                    results = await asyncio.gather(*tasks, return_exceptions=True)
                    all_results.extend([r for r in results 
                                      if isinstance(r, AcquisitionResult)])
                
                return rcsb_id, all_results

        # Process all structures
        tasks = [process_structure(rcsb_id) for rcsb_id in rcsb_ids]
        results = await asyncio.gather(*tasks, return_exceptions=True)
        
        return {
            rcsb_id: results 
            for rcsb_id, results in results 
            if isinstance(results, tuple)
        }

    return asyncio.run(_process_chunk_async())