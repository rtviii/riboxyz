from functools import partial
import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.global_ops import GlobalOps
from ribctl.ribosome_ops import RibosomeOps
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from ribctl import RIBETL_DATA
from ribctl.asset_manager.parallel_acquisition import process_chunk
import click
from typing import List, Tuple
import multiprocessing
from concurrent.futures import ALL_COMPLETED, Future, ProcessPoolExecutor, ThreadPoolExecutor, wait
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

@db.command()
@click.argument('pdb_ids', type=PDBIDsParam(), nargs=-1)
@click.option('--workers', '-w',
              default=10,
              help='Number of worker threads')
@click.option('--force', '-f', is_flag=True, help='Force upload even if structure exists')
@click.option('--mode', '-m', type=click.Choice(['full', 'structure', 'ligands'], case_sensitive=False),
               default='full',
              help='Upload mode: full (all data), structure (only structure nodes), or ligands')
@click.option('--instance', '-i', help='Neo4j database instance name', default=NEO4J_CURRENTDB)
@click.pass_context
def upload(ctx, pdb_ids, workers, force, mode, instance):
    """Upload structure data to Neo4j database.
    
    Examples:\n
    \b
    # Upload complete data for specific structures
    ribd.py db upload 3J7Z 4V6X
    
    \b
    # Upload only structure nodes to specific instance
    ribd.py db upload --mode structure --instance neo4j 3J7Z
    
    \b
    # Upload from file with force option
    cat pdb_list.txt | ribd.py db upload --force
    """
    all_pdb_ids = list(set(ctx.obj['piped_pdb_ids'] + sum(pdb_ids, [])))
    
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return
    
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance, NEO4J_PASSWORD)
    except Exception as e:
        click.echo(f"Failed to connect to database: {str(e)}", err=True)
        return

    click.echo(f"Starting upload to instance '{instance}'...")
    with click.progressbar(length=len(all_pdb_ids),
                         label='Uploading structures',
                         show_pos=True,
                         file=sys.stderr) as bar:
        
        futures: list[Future] = []
        
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for rcsb_id in sorted(all_pdb_ids):
                if mode == 'full':
                    fut = executor.submit(
                        partial(adapter.add_total_structure, rcsb_id, force)
                    )
                elif mode == 'structure':
                    fut = executor.submit(
                        partial(adapter.upsert_structure_node, rcsb_id)
                    )
                elif mode == 'ligands':
                    profile = RibosomeOps(rcsb_id).profile
                    for ligand in profile.nonpolymeric_ligands:
                        if not "ion" in ligand.chemicalName.lower():
                            fut = executor.submit(
                                partial(adapter.upsert_ligand_node, ligand, rcsb_id)
                            )
                
                fut.add_done_callback(lambda p: bar.update(1))
                futures.append(fut)
                
            # Wait for all uploads to complete
            errors = []
            done, _ = wait(futures, return_when=ALL_COMPLETED)
            for future in done:
                try:
                    future.result()  # Check for exceptions
                except Exception as e:
                    errors.append(str(e))

    if errors:
        click.echo("\nErrors occurred during upload:", err=True)
        for error in errors:
            click.echo(f"  - {error}", err=True)
    else:
        click.echo("\nUpload completed successfully")

@db.command()
@click.option('--force', '-f', is_flag=True,
              help='Force reinitialization of database')
@click.argument('instance_name', required=True)
@click.pass_context
def init(ctx, force, instance_name):
    """Initialize a new Neo4j database instance.
    
    Creates necessary constraints and initial data structures in the specified instance.
    
    Examples:\n
    \b
    # Initialize a new instance named 'ribosome'
    ribd.py db init ribosome
    
    \b
    # Force reinitialize an existing instance
    ribd.py db init --force ribosome
    """
    if not force:
        click.confirm(f'This will initialize database instance "{instance_name}". Continue?',
                     abort=True)
    
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance_name, NEO4J_PASSWORD)
        adapter.initialize_new_instance()
        click.echo(f"Database instance '{instance_name}' initialized successfully")
    except Exception as e:
        click.echo(f"Failed to initialize database instance: {str(e)}", err=True)

@db.group()
@click.pass_context
def instance(ctx):
    """Manage Neo4j database instances"""
    pass

@instance.command(name='create')
@click.argument('name', required=True)
@click.option('--force', '-f', is_flag=True,
              help='Force creation if instance exists')
@click.pass_context
def create_instance(ctx, name, force):
    """Create a new Neo4j database instance.
    
    Examples:\n
    \b
    # Create a new instance
    ribd.py db instance create ribosome2
    
    \b
    # Force create an instance
    ribd.py db instance create --force ribosome2
    """
    try:
        # Create database instance using system database
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, "system", NEO4J_PASSWORD)
        with adapter.driver.session() as session:
            # Check if instance exists
            result = session.run("SHOW DATABASES")
            databases = [record["name"] for record in result]
            
            if name in databases and not force:
                click.echo(f"Database instance '{name}' already exists. Use --force to recreate.", err=True)
                return
            
            if force:
                session.run(f"DROP DATABASE {name} IF EXISTS")
            
            # Create new database
            session.run(f"CREATE DATABASE {name}")
            click.echo(f"Created database instance '{name}'")
            
            # Initialize the new instance
            adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, name, NEO4J_PASSWORD)
            adapter.initialize_new_instance()
            click.echo(f"Initialized database instance '{name}'")
            
    except Exception as e:
        click.echo(f"Failed to create database instance: {str(e)}", err=True)


@db.command(name='upload_all')
@click.option('--workers', '-w',
              default=10,
              help='Number of worker threads')
@click.option('--force', '-f', is_flag=True, 
              help='Force upload even if structures exist')
@click.option('--mode', '-m', 
              type=click.Choice(['full', 'structure', 'ligands'], case_sensitive=False),
              default='full',
              help='Upload mode: full (all data), structure (only structure nodes), or ligands')
@click.option('--instance', '-i', 
              help='Neo4j database instance name', 
              default=NEO4J_CURRENTDB)
@click.option('--chunk-size', '-c', 
              default=10,
              help='Number of structures to process per chunk')
@click.pass_context
def upload_all(ctx, workers, force, mode, instance, chunk_size):
    """Upload all available structures to Neo4j database in parallel.
    
    Examples:\n
    \b
    # Upload all structures with default settings
    ribd.py db upload_all
    
    \b
    # Upload all structures to specific instance with more workers
    ribd.py db upload_all --instance neo4j --workers 20
    
    \b
    # Upload only structure nodes for all structures
    ribd.py db upload_all --mode structure
    """
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance, NEO4J_PASSWORD)
    except Exception as e:
        click.echo(f"Failed to connect to database: {str(e)}", err=True)
        return

    # Get all available structure IDs
    try:
        all_structures = GlobalOps.list_all_structs()
        if not all_structures:
            click.echo("No structures found to upload", err=True)
            return
        
        click.echo(f"Found {len(all_structures)} structures to process")
    except Exception as e:
        click.echo(f"Failed to get structure list: {str(e)}", err=True)
        return

    # Split structures into chunks for parallel processing
    chunks = [all_structures[i:i + chunk_size] 
             for i in range(0, len(all_structures), chunk_size)]

    click.echo(f"Starting parallel upload to instance '{instance}'...")
    click.echo(f"Processing in chunks of {chunk_size} with {workers} workers")

    with tqdm(total=len(all_structures), 
             desc="Uploading structures",
             unit="structure") as pbar:
        
        futures: list[Future] = []
        errors: list[tuple[str, str]] = []  # (structure_id, error_message)
        
        with ThreadPoolExecutor(max_workers=workers) as executor:
            for chunk in chunks:
                for rcsb_id in chunk:
                    if mode == 'full':
                        fut = executor.submit(
                            partial(adapter.add_total_structure, rcsb_id, force)
                        )
                    elif mode == 'structure':
                        fut = executor.submit(
                            partial(adapter.upsert_structure_node, rcsb_id)
                        )
                    elif mode == 'ligands':
                        try:
                            profile = RibosomeOps(rcsb_id).profile
                            for ligand in profile.nonpolymeric_ligands:
                                if not "ion" in ligand.chemicalName.lower():
                                    fut = executor.submit(
                                        partial(adapter.upsert_ligand_node, ligand, rcsb_id)
                                    )
                        except Exception as e:
                            errors.append((rcsb_id, f"Failed to process ligands: {str(e)}"))
                            continue

                    def update_progress(future, rcsb_id=rcsb_id):
                        try:
                            future.result()  # Check for exceptions
                        except Exception as e:
                            errors.append((rcsb_id, str(e)))
                        pbar.update(1)
                        pbar.set_description(f"Processed {rcsb_id}")

                    fut.add_done_callback(update_progress)
                    futures.append(fut)

            # Wait for all uploads to complete
            done, _ = wait(futures, return_when=ALL_COMPLETED)

    # Report results
    total_processed = len(all_structures)
    failed = len(errors)
    succeeded = total_processed - failed

    click.echo("\nUpload Summary:")
    click.echo(f"Total structures processed: {total_processed}")
    click.echo(f"Successfully uploaded: {succeeded}")
    click.echo(f"Failed: {failed}")

    if errors:
        click.echo("\nErrors occurred during upload:")
        for rcsb_id, error in errors:
            click.echo(f"  - {rcsb_id}: {error}")
    else:
        click.echo("\nAll structures uploaded successfully")


















if __name__ == '__main__':
    cli(obj={})
