import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl import RIBETL_DATA
from ribctl.asset_manager.parallel_acquisition import process_chunk
import click
from typing import List, Tuple
from pathlib import Path
import multiprocessing
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm

from ribctl.asset_manager.asset_types import AssetType

def get_input_pdb_ids() -> List[str]:
    """Get PDB IDs from either stdin (if piped) or return None to handle as argument"""
    if not sys.stdin.isatty():
        return [line.strip().upper() for line in sys.stdin if line.strip()]
    return []

class PDBIDsParam(click.ParamType):
    """Custom parameter type for PDB IDs validation"""
    name = "pdb_ids"
    
    def convert(self, value, param, ctx):
        if isinstance(value, list):
            return value
        if not value:
            return []
        pdb_id = value.strip().upper()
        if not (len(pdb_id) == 4 and pdb_id.isalnum()):
            self.fail(f"Invalid PDB ID format: {value}", param, ctx)
        return [pdb_id]

@click.group()
@click.pass_context
def cli(ctx):
    """ribctl - Command line interface for ribosome data pipeline"""
    ctx.ensure_object(dict)
    ctx.obj['piped_pdb_ids'] = get_input_pdb_ids()

@cli.group()
@click.pass_context
def etl(ctx):
    """ETL operations for ribosome data"""
    pass

@cli.group()
@click.pass_context
def db(ctx):
    """Database operations for ribosome data"""
    pass

@etl.command()
@click.argument('pdb_ids', type=PDBIDsParam(), nargs=-1)
@click.pass_context
def verify(ctx, pdb_ids):
    """Verify assets exist for given PDB IDs"""
    all_pdb_ids = list(set(ctx.obj['piped_pdb_ids'] + sum(pdb_ids, [])))
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return
    
    click.echo(f"Verifying assets for: {', '.join(all_pdb_ids)}")
    # TODO: Implement verification logic

@etl.command()
@click.pass_context
def verify_all(ctx):
    """Verify all assets in the system"""
    click.echo("Verifying all assets in the system...")
    # TODO: Implement full verification logic
    # This should probably:
    # 1. Get list of all expected assets
    # 2. Check their existence and integrity
    # 3. Report any issues found

@etl.command()
@click.option('--asset-types', '-t', 
              multiple=True,
              type=click.Choice([t.name for t in AssetType], case_sensitive=True),
              help='Asset types to acquire (if not specified, all assets will be acquired)')
@click.option('--force', '-f', is_flag=True,
              help='Force regeneration of existing assets')
@click.option('--workers', '-w',
              default=max(1, multiprocessing.cpu_count() - 1),
              help='Number of worker processes')
@click.option('--chunk-size', '-c', default=4,
              help='Number of structures to process per worker')
@click.option('--concurrent-structures', '-s', default=4,
              help='Maximum concurrent structures per worker')
@click.option('--concurrent-assets', '-a', default=3,
              help='Maximum concurrent assets per structure')
@click.argument('pdb_ids', type=PDBIDsParam(), nargs=-1)
@click.pass_context
def get(ctx, asset_types, force, workers, chunk_size,
        concurrent_structures, concurrent_assets, pdb_ids):
    """Download or generate assets for given PDB IDs. 
    
    Examples:\n
    \b
    # Get specific assets for a single structure
    ribd.py etl get --asset-types MMCIF STRUCTURE_PROFILE 3J7Z
    
    \b
    # Get all assets for multiple structures
    ribd.py etl get 3J7Z 4V6X
    
    \b
    # Process structures from a file with specific assets
    cat pdb_list.txt | ribd.py etl get --asset-types MMCIF PTC
    """
    all_pdb_ids = list(set(ctx.obj['piped_pdb_ids'] + sum(pdb_ids, [])))
    
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return

    # If no asset types specified, use all available types
    selected_asset_types = (
        [AssetType[name] for name in asset_types]
        if asset_types
        else list(AssetType)
    )

    # Inform user about the operation
    assets_str = ", ".join(ast.name for ast in selected_asset_types)
    click.echo(f"Getting assets [{assets_str}] for {len(all_pdb_ids)} structures...")

    # Split PDB IDs into chunks
    chunks = [all_pdb_ids[i:i + chunk_size] 
             for i in range(0, len(all_pdb_ids), chunk_size)]

    with tqdm(total=len(all_pdb_ids)) as pbar:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = [
                executor.submit(
                    process_chunk,
                    str(RIBETL_DATA),
                    chunk,
                    [ast.name for ast in selected_asset_types],  # Pass enum names
                    force,
                    concurrent_structures,
                    concurrent_assets
                )
                for chunk in chunks
            ]

            # Process results as they complete
            for future in futures:
                try:
                    results = future.result()
                    for rcsb_id, asset_results in results.items():
                        pbar.update(1)
                        pbar.set_description(f"Processed {rcsb_id}")
                        
                        # Report failures
                        for result in asset_results:
                            if not result.success:
                                click.echo(
                                    f"Error with {result.asset_type_name} "
                                    f"for {rcsb_id}: {result.error}",
                                    err=True
                                )
                except Exception as e:
                    click.echo(f"Error processing chunk: {str(e)}", err=True)

@db.command()
@click.argument('pdb_ids', type=PDBIDsParam(), nargs=-1)
@click.option('--asset-type', '-t', 
              type=click.Choice(['all', 'structure', 'sequence', 'alignment']),
              default='all', 
              help='Type of assets to upload')
@click.pass_context
def upload(ctx, pdb_ids, asset_type):
    """Upload assets to database for given PDB IDs"""
    all_pdb_ids = list(set(ctx.obj['piped_pdb_ids'] + sum(pdb_ids, [])))
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return
        
    click.echo(f"Uploading {asset_type} assets for: {', '.join(all_pdb_ids)}")
    # TODO: Implement upload logic

@db.command()
@click.argument('query_file', type=click.Path(exists=True))
@click.option('--params', '-p', multiple=True, help='Parameters for the Cypher query in key=value format')
@click.option('--output', '-o', type=click.Path(), help='Output file for query results')
@click.pass_context
def cypher(ctx, query_file, params, output):
    """Execute a Cypher query from a file"""
    # Parse parameters
    query_params = {}
    for param in params:
        try:
            key, value = param.split('=')
            query_params[key.strip()] = value.strip()
        except ValueError:
            click.echo(f"Invalid parameter format: {param}. Use key=value format.", err=True)
            return

    # Read query from file
    with open(query_file, 'r') as f:
        query = f.read()

    click.echo(f"Executing Cypher query from {query_file}")
    # TODO: Implement query execution logic
    # results = execute_cypher_query(query, query_params)

    # if output:
    #     # Save results to file
    #     with open(output, 'w') as f:
    #         json.dump(results, f)
    # else:
    #     # Print results to stdout
    #     click.echo(results)

if __name__ == '__main__':
    cli(obj={})