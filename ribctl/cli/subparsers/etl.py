import asyncio
from enum import Enum
from typing import List, Set
import click
from click import Context
from loguru import logger
from pathlib import Path

from ribctl import RIBETL_DATA
from ribctl.etl.asset_manager import AssetType, RibosomeAssetManager
from ribctl.etl.assets_global import GlobalAssets
from ribctl.logs.loggers import get_etl_logger

# Initialize the asset manager globally
ASSET_MANAGER = GlobalAssets(RIBETL_DATA)

def format_asset_table(status_dict: dict) -> List[str]:
    """Format asset status into a pretty table"""
    # Header
    headers = ["RCSB_ID"] + [asset.name for asset in AssetType]
    header_row = "\t".join(headers)
    
    # Data rows
    rows = [header_row]
    for pdb_id, status in status_dict.items():
        row = [pdb_id] + ['X' if status.get(asset, False) else ' ' for asset in AssetType]
        rows.append("\t".join(row))
    
    return rows

@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    """Asset management CLI for ribosome structures"""
    ctx.ensure_object(dict)

@etl.command()
@click.pass_context
def status(ctx: Context):
    """Display status of all assets for all structures"""
    async def get_all_status():
        structures = ASSET_MANAGER.list_all_structs()
        all_status = {}
        
        for pdb_id in structures:
            status = await ASSET_MANAGER.verify_all_assets(pdb_id)
            all_status[pdb_id] = status
            
        return all_status
    
    # Get and display status
    status_dict = asyncio.run(get_all_status())
    rows = format_asset_table(status_dict)
    
    # Print with headers every 20 rows
    for i, row in enumerate(rows):
        if i == 0 or i % 20 == 0:
            click.echo(rows[0])  # Print header
        click.echo(row)

@etl.command()
@click.pass_context
@click.argument(
    "assets",
    required=True,
    nargs=-1,
    type=click.Choice([t.name for t in AssetType])
)
@click.option("--overwrite", is_flag=True, default=False, help="Overwrite existing assets")
@click.option("--reclassify", is_flag=True, default=False, help="Reclassify structure")
@click.option("--rcsb-sync", is_flag=True, default=False, help="Sync with RCSB")
@click.option("--all-structs", is_flag=True, default=False, help="Process all structures")
def generate(
    ctx: Context,
    assets: tuple[str],
    overwrite: bool,
    reclassify: bool,
    rcsb_sync: bool,
    all_structs: bool
):
    """Generate or update assets for structures"""
    
    # Convert string asset names to AssetType enum members
    asset_types = [AssetType[asset] for asset in assets]
    
    # async def process_structure(pdb_id: str, asset_types: List[AssetType]):
    #     """Process assets for a single structure"""
    #     try:
    #         for asset_type in asset_types:
    #             await ASSET_MANAGER.generate_asset(
    #                 pdb_id,
    #                 asset_type,
    #                 force=overwrite
    #             )
    #         logger.info(f"Processed successfully {pdb_id}: {[a.name for a in asset_types]}")
    #     except Exception as e:
    #         logger.error(f"Error processing {pdb_id}: {e}")
    #         click.echo(f"Error processing {pdb_id}: {e}", err=True)

    async def main():
        if 'rcsb_id' in ctx.obj:
            # Single structure processing
            await process_structure(ctx.obj['rcsb_id'].upper(), asset_types)
            return

        if rcsb_sync:
            # Sync with RCSBwj
            structures =  ASSET_MANAGER.status_vs_rcsb()
            sync_assets = [AssetType.MMCIF, AssetType.STRUCTURE_PROFILE]
            for pdb_id in structures:
                click.echo(f"RCSB Sync: Fetching assets for {pdb_id}")
                # await process_structure(pdb_id, sync_assets)
            return

        if all_structs:
            # Process all structures
            structures =  ASSET_MANAGER.list_all_structs()
            for pdb_id in structures:
                click.echo(f"[{pdb_id}]")
                # await process_structure(pdb_id, asset_types)
            return

    asyncio.run(main())

@etl.command()
@click.pass_context
@click.argument("pdb_id", required=True)
@click.argument(
    "asset",
    required=True,
    type=click.Choice([t.name for t in AssetType])
)
def verify(ctx: Context, pdb_id: str, asset: str):
    """Verify specific asset for a structure"""
    async def verify_asset():
        asset_type = AssetType[asset]
        is_valid = await ASSET_MANAGER.verify_asset(pdb_id, asset_type)
        status = "✓" if is_valid else "✗"
        click.echo(f"{pdb_id} {asset}: {status}")
    
    asyncio.run(verify_asset())

if __name__ == "__main__":
    etl(obj={})