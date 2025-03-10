import os
from pprint import pprint
import sys
from loguru import logger

sys.dont_write_bytecode = True
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.lib.libbsite import extract_ligand_to_mmcif
from logger_config import configure_logging

configure_logging()
import concurrent.futures as cf
from neo4j_ribosome.db_lib_reader import Neo4jReader
from functools import partial
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_registry import AssetRegistry
from ribctl.lib.schema.types_ribosome import PTCInfo, RibosomeStructure
from ribctl.asset_manager.assets_structure import StructureAssets
from ribctl.global_ops import GlobalOps
from ribctl.ribosome_ops import RibosomeOps
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from ribctl import RIBETL_DATA
from typing import List, Optional, Tuple, Dict
import click
from pathlib import Path
import asyncio
from concurrent.futures import ProcessPoolExecutor
import math
from ribctl.asset_manager.parallel_acquisition import (
    process_chunk,
    AcquisitionResult,
    process_chunk_with_tracking,
)
import multiprocessing
from ribctl.asset_manager.asset_registry import main_registry
from concurrent.futures import (
    ProcessPoolExecutor,
)
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

        if type(value) == list:  # Changed from isinstance(value, list)
            return value
        if not value:
            return []
        pdb_id = value.strip().upper()
        if not (len(pdb_id) == 4 and pdb_id.isalnum()):
            self.fail(f"Invalid PDB ID format: {value}", param, ctx)
        return [pdb_id]


def get_input_structures() -> List[str]:
    """Read PDB IDs from stdin if available, otherwise return empty list"""
    if not sys.stdin.isatty():
        return [line.strip().upper() for line in sys.stdin.readlines() if line.strip()]
    return []


@click.group()
@click.pass_context
def cli(ctx):
    """ribctl - Command line interface for ribosome data pipeline"""
    ctx.ensure_object(dict)
    ctx.obj["piped_pdb_ids"] = get_input_pdb_ids()


@cli.group()
@click.pass_context
def etl(ctx):
    """ETL operations for ribosome data"""
    pass


@etl.command()
@click.pass_context
@click.option(
    "-t",
    "--asset-type",
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    multiple=True,
    required=True,
    help="Asset type(s) to acquire. Can be specified multiple times.",
)
@click.argument("pdb_id", required=False)
@click.option("--force", is_flag=True, help="Force regeneration of existing assets")
@click.option(
    "--max-structures",
    default=4,
    help="Maximum number of concurrent structures to process",
)
@click.option(
    "--max-assets", default=3, help="Maximum number of concurrent assets per structure"
)
def get(
    ctx,
    asset_type: List[str],
    pdb_id: str,
    force: bool,
    max_structures: int,
    max_assets: int,
):
    """
    Acquire assets for specified structure(s). Accepts either a single PDB ID as argument
    or a list of IDs via stdin pipe.
    """
    # Get structures from either argument or pipe
    structures = ctx.obj.get("piped_pdb_ids", [])
    
    # If nothing in context, try reading from stdin directly
    if not structures:
        structures = get_input_structures()
        
    if pdb_id:
        structures.append(pdb_id.upper())
    

    # Convert asset type names to enum
    asset_types = [AssetType[t] for t in asset_type]

    async def process_structure(rcsb_id: str):
        try:
            await main_registry.generate_multiple(rcsb_id, asset_types, force)
            click.echo(f"Successfully processed assets for {rcsb_id}")
        except Exception as e:
            logger.exception(f"Failed to process {rcsb_id}")
            click.echo(f"Failed to process {rcsb_id}: {str(e)}", err=True)

    async def process_all():
        sem = asyncio.Semaphore(max_structures)

        async def wrapped_process(rcsb_id: str):
            async with sem:
                await process_structure(rcsb_id)

        tasks = [wrapped_process(rcsb_id) for rcsb_id in structures]
        await asyncio.gather(*tasks)

    # Run the async processing
    asyncio.run(process_all())


@etl.command()
@click.option(
    "-t",
    "--asset-type",
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    multiple=True,
    required=True,
    help="Asset type(s) to acquire. Can be specified multiple times.",
)
@click.option("--force", is_flag=True, help="Force regeneration of existing assets")
@click.option(
    "--max-structures",
    default=4,
    help="Maximum number of concurrent structures per process",
)
@click.option(
    "--max-assets", default=3, help="Maximum number of concurrent assets per structure"
)
@click.option("--processes", default=4, help="Number of parallel processes to use")
@click.option(
    "--chunk-size", default=10, help="Number of structures to process per chunk"
)
def get_all(
    asset_type: List[str],
    force: bool,
    max_structures: int,
    max_assets: int,
    processes: int,
    chunk_size: int,
):
    """
    Acquire assets for ALL structures in the database using parallel processing.
    Structures are processed in parallel chunks using multiple CPU cores.
    """
    # Convert asset type names to enum
    asset_types = [AssetType[t] for t in asset_type]

    # Get all available structures
    structures = GlobalOps.list_profiles()

    if not structures:
        click.echo("No structures found in the database!", err=True)
        return

    total_structures = len(structures)
    click.echo(f"Found {total_structures} structures to process")

    # Split structures into chunks
    chunks = [
        structures[i : i + chunk_size] for i in range(0, len(structures), chunk_size)
    ]

    # Create a partial function with all the fixed arguments
    process_func = partial(
        process_chunk_with_tracking,
        base_dir=str(Path(RIBETL_DATA)),
        asset_type_names=[t.name for t in asset_types],
        force=force,
        max_structures=max_structures,
        max_assets=max_assets,
    )

    # Process chunks in parallel
    click.echo(f"Processing {len(chunks)} chunks using {processes} processes")
    processed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=processes) as executor:
        for result in executor.map(process_func, chunks):
            # Count successes and failures from results
            for rcsb_id, acquisitions in result.items():
                if all(acq.success for acq in acquisitions):
                    processed += 1
                else:
                    failed += 1

    # Print final summary
    click.echo(f"\nProcessing complete:")
    click.echo(f"Total structures: {total_structures}")
    click.echo(f"Successfully processed: {processed}")
    click.echo(f"Failed: {failed}")

    # Return non-zero exit code if any failures
    if failed > 0:
        raise click.ClickException(f"Failed to process {failed} structures")


@etl.command()
@click.argument("pdb_ids", type=PDBIDsParam(), nargs=-1)
@click.pass_context
def verify(ctx, pdb_ids):
    """Verify assets exist for given PDB IDs"""
    all_pdb_ids = list(set(ctx.obj["piped_pdb_ids"] + sum(pdb_ids, [])))
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return




@etl.command()
@click.pass_context
def verify_all(ctx):
    """Verify all assets in the system by validating each structure against the RibosomeStructure model"""
    from concurrent.futures import ThreadPoolExecutor
    import asyncio
    from typing import Dict, List, Tuple

    structures = GlobalOps.list_profiles()
    if not structures:
        click.echo("No structures found in the system", err=True)
        return

    click.echo(f"Found {len(structures)} structures to verify")

    # Track validation results
    results: Dict[str, List[Tuple[str, str]]] = {
        "valid": [],
        "invalid": []
    }

    async def validate_structure(rcsb_id: str) -> None:
        """Validate a single structure and update results"""
        try:
            ops = RibosomeOps(rcsb_id)
            structure_data = ops.assets.profile()
            
            # Attempt validation
            RibosomeStructure.model_validate(structure_data)
            results["valid"].append(rcsb_id)
            
        except Exception as e:
            results["invalid"].append((rcsb_id, str(e)))

    async def process_all_structures():
        """Process all structures concurrently with a progress bar"""
        tasks = []
        with click.progressbar(
            structures,
            label="Validating structures",
            length=len(structures)
        ) as progress_structures:
            for rcsb_id in progress_structures:
                task = asyncio.create_task(validate_structure(rcsb_id))
                tasks.append(task)
            await asyncio.gather(*tasks)

    # Run validation
    try:
        asyncio.run(process_all_structures())
    except Exception as e:
        click.echo(f"\nError during validation: {str(e)}", err=True)
        return

    # Report results
    total = len(structures)
    valid_count = len(results["valid"])
    invalid_count = len(results["invalid"])

    click.echo("\nValidation Summary:")
    click.echo(f"Total structures processed: {total}")
    click.echo(f"Valid structures: {valid_count}")
    click.echo(f"Invalid structures: {invalid_count}")

    if invalid_count > 0:
        click.echo("\nInvalid structures:")
        for rcsb_id, error in results["invalid"]:
            click.echo(f"  - {rcsb_id}: {error}")

    # Return non-zero exit code if any structures failed validation
    if invalid_count > 0:
        ctx.exit(1)



@etl.command()
@click.pass_context
def ligand_mmcifs(ctx):
    from neo4j_ribosome.db_lib_reader import dbqueries

    ligs = dbqueries.ligands_per_structure()
    for lig in ligs:
        chemid = lig["chemicalId"]
        structs = lig["structures"]
        for struct in structs:
            structure = RibosomeOps(struct).assets.biopython_structure()
            try:
                if os.path.exists(StructureAssets(struct).paths.nonpoly_entity(chemid)):
                    print("Already exists: ", struct, chemid)
                    continue
                extract_ligand_to_mmcif(
                    structure,
                    chemid,
                    StructureAssets(struct).paths.nonpoly_entity(chemid),
                )
            except:
                pass
# Define this outside the command function to make it picklable
# For multiprocessing
import concurrent.futures
from tqdm import tqdm
import sys
import time

def validate_structure(rcsb_id):
    """Validate a single structure against the Pydantic model"""
    from pydantic import ValidationError
    
    try:
        ops = RibosomeOps(rcsb_id)
        # Get structure profile
        structure_data = ops.assets.profile()
        
        # Validate against Pydantic model
        RibosomeStructure.model_validate(structure_data)
        
        return {
            "rcsb_id": rcsb_id,
            "valid": True,
            "error": None
        }
    except ValidationError as e:
        return {
            "rcsb_id": rcsb_id,
            "valid": False,
            "error": str(e)
        }
    except Exception as e:
        return {
            "rcsb_id": rcsb_id,
            "valid": False,
            "error": f"Failed to load structure: {str(e)}"
        }


@etl.command()
@click.option(
    "--show-errors",
    is_flag=True,
    help="Show detailed Pydantic validation errors for failed structures",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Show detailed progress information",
)
@click.option(
    "--workers",
    default=max(1, multiprocessing.cpu_count() - 1),
    help="Number of worker processes for parallel validation",
)
def check_schema(show_errors, verbose, workers):
    """
    Validate all structure profiles against the Pydantic model schema.
    
    Attempts to load each structure profile and validates it against the
    RibosomeStructure Pydantic model. Reports structures that fail validation.
    
    Examples:\n
    \b
    # Check all structures with minimal output
    ribd.py etl check_schema
    
    \b
    # Show detailed validation errors
    ribd.py etl check_schema --show-errors
    
    \b
    # Verbose output with 8 worker processes
    ribd.py etl check_schema -v --workers 8
    """
    from concurrent.futures import ProcessPoolExecutor
    import time
    import concurrent.futures
    
    # Get all available structures
    structures = GlobalOps.list_profiles()
    total = len(structures)
    
    if not structures:
        click.echo("No structures found.")
        return
    
    click.echo(f"Validating schema for {total} structures...")
    
    # Results tracking
    results = {
        "passed": [],
        "failed": [],
        "errors": {}
    }
    
    start_time = time.time()
    
    # Process structures in parallel
    with ProcessPoolExecutor(max_workers=workers) as executor:
        with tqdm(total=total, disable=not verbose) as progress_bar:
            futures = [executor.submit(validate_structure, rcsb_id) for rcsb_id in structures]
            
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                rcsb_id = result["rcsb_id"]
                
                if result["valid"]:
                    results["passed"].append(rcsb_id)
                    if verbose:
                        click.echo(f"✓ {rcsb_id}: Valid")
                else:
                    results["failed"].append(rcsb_id)
                    results["errors"][rcsb_id] = result["error"]
                    if verbose:
                        click.echo(f"✗ {rcsb_id}: Invalid")
                
                progress_bar.update(1)
    
    # Print summary
    elapsed_time = time.time() - start_time
    passed_count = len(results["passed"])
    failed_count = len(results["failed"])
    
    click.echo("\nValidation Summary:")
    click.echo(f"Total structures: {total}")
    click.echo(f"Passed: {passed_count} ({passed_count/total*100:.1f}%)")
    click.echo(f"Failed: {failed_count} ({failed_count/total*100:.1f}%)")
    click.echo(f"Time elapsed: {elapsed_time:.2f} seconds")
    
    # Print failed structures
    if failed_count > 0:
        click.echo("\nFailed structures:")
        for i, rcsb_id in enumerate(sorted(results["failed"]), 1):
            click.echo(f"{rcsb_id}")
            
            # Print detailed errors if requested
            if show_errors:
                error_msg = results["errors"][rcsb_id]
                formatted_error = "\n".join(f"    {line}" for line in error_msg.split("\n"))
                click.echo(f"   Error details:\n{formatted_error}\n")
    
    # Return non-zero exit code if any structures failed validation
    if failed_count > 0:
        sys.exit(1)

@etl.command(name="sync_all")
@click.option(
    "--asset-types",
    "-t",
    multiple=True,
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    help="Asset types to acquire (required)",
    required=True,
)
@click.option(
    "--force", "-f", is_flag=True, help="Force regeneration of existing assets"
)
@click.option(
    "--workers",
    "-w",
    default=max(1, multiprocessing.cpu_count() - 1),
    help="Number of worker processes",
)
@click.option(
    "--chunk-size", "-c", default=4, help="Number of structures to process per worker"
)
@click.option(
    "--concurrent-structures",
    "-s",
    default=4,
    help="Maximum concurrent structures per worker",
)
@click.option(
    "--concurrent-assets",
    "-a",
    default=3,
    help="Maximum concurrent assets per structure",
)
@click.option("--delay", "-d", default=2.0, help="Delay in seconds between chunks")
@click.option(
    "--max-retries",
    "-r",
    default=3,
    help="Maximum number of retries per failed operation",
)
@click.pass_context
def sync_all(
    ctx,
    asset_types,
    force,
    workers,
    chunk_size,
    concurrent_structures,
    concurrent_assets,
    delay,
    max_retries,
):
    """Synchronize all structure assets with RCSB, downloading or updating as needed.

    Compares local assets with RCSB database and only processes structures that:
    - Don't exist locally
    - Have different modification dates
    - Are missing requested assets

    Examples:\n
    \b
    # Sync specific assets for all outdated structures
    ribd.py etl sync_all --asset-types MMCIF STRUCTURE_PROFILE

    \b
    # Force sync with custom settings
    ribd.py etl sync_all --asset-types MMCIF PTC --workers 4 --force
    """
    from time import sleep
    import random

    def process_with_retry(operation, max_retries):
        """Execute operation with retries and exponential backoff"""
        for attempt in range(max_retries):
            try:
                return operation()
            except Exception as e:
                if attempt == max_retries - 1:
                    raise e
                sleep_time = (attempt + 1) * 2 + random.random()
                sleep(sleep_time)

    try:
        # Get structures that need updating
        click.echo("Comparing local assets with RCSB database...")

        # Here we'd use your vs_rcsb method to get structures needing updates
        structures_to_update = GlobalOps.missing_profiles()

        if not structures_to_update:
            click.echo("All structures are up to date!")
            return

        click.echo(f"Found {len(structures_to_update)} structures needing updates")
    except Exception as e:
        click.echo(f"Failed to compare with RCSB: {str(e)}", err=True)
        return

    # Convert asset type names to enums
    selected_asset_types = [AssetType[name] for name in asset_types]
    assets_str = ", ".join(ast.name for ast in selected_asset_types)

    # Inform user about the operation
    click.echo(
        f"Getting assets [{assets_str}] for {len(structures_to_update)} structures..."
    )
    click.echo(f"Using {workers} workers, {chunk_size} structures per chunk")
    click.echo(f"Delay between chunks: {delay}s")

    # Split structures into chunks
    chunks = [
        structures_to_update[i : i + chunk_size]
        for i in range(0, len(structures_to_update), chunk_size)
    ]

    errors = []
    with tqdm(total=len(structures_to_update)) as pbar:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = []

            # Submit all chunks to the process pool
            for chunk in chunks:
                future = executor.submit(
                    process_chunk,
                    str(RIBETL_DATA),
                    chunk,
                    [ast.name for ast in selected_asset_types],
                    force,
                    concurrent_structures,
                    concurrent_assets,
                )
                futures.append(future)

            # Process results as they complete
            for idx, future in enumerate(futures):
                try:
                    results = process_with_retry(lambda: future.result(), max_retries)

                    # Process results for this chunk
                    for rcsb_id, asset_results in results.items():
                        pbar.update(1)
                        pbar.set_description(f"Processed {rcsb_id}")

                        # Report failures
                        for result in asset_results:
                            if not result.success:
                                error_msg = (
                                    f"Error with {result.asset_type_name} "
                                    f"for {rcsb_id}: {result.error}"
                                )
                                errors.append((rcsb_id, error_msg))
                                click.echo(f"\n{error_msg}", err=True)

                    # Apply delay after each chunk except the last
                    if idx < len(futures) - 1:
                        sleep(delay)

                except Exception as e:
                    error_msg = (
                        f"Failed to process chunk after {max_retries} retries: {str(e)}"
                    )
                    click.echo(f"\n{error_msg}", err=True)
                    chunk = chunks[idx]  # Get corresponding chunk for this future
                    errors.extend((str(pdb_id), error_msg) for pdb_id in chunk)

    # Final report
    total_processed = len(structures_to_update)
    failed = len(errors)
    succeeded = total_processed - failed

    click.echo("\nSync Summary:")
    click.echo(f"Total structures processed: {total_processed}")
    click.echo(f"Successfully processed: {succeeded}")
    click.echo(f"Failed: {failed}")

    if errors:
        click.echo("\nErrors occurred during processing:")
        for rcsb_id, error in errors:
            click.echo(f"  - {rcsb_id}: {error}")
    else:
        click.echo("\nAll structures processed successfully")


@cli.group()
@click.pass_context
def db(ctx):
    """Database operations for ribosome data"""
    pass


@db.command()
@click.argument("query_file", type=click.Path(exists=True))
@click.option(
    "--params",
    "-p",
    multiple=True,
    help="Parameters for the Cypher query in key=value format",
)
@click.option("--output", "-o", type=click.Path(), help="Output file for query results")
@click.pass_context
def cypher(ctx, query_file, params, output):
    """Execute a Cypher query from a file"""
    # Parse parameters
    query_params = {}
    for param in params:
        try:
            key, value = param.split("=")
            query_params[key.strip()] = value.strip()
        except ValueError:
            click.echo(
                f"Invalid parameter format: {param}. Use key=value format.", err=True
            )
            return

    # Read query from file
    with open(query_file, "r") as f:
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
@click.option("--force", "-f", is_flag=True, help="Force reinitialization of database")
@click.argument("instance_name", required=True)
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
        click.confirm(
            f'This will initialize database instance "{instance_name}". Continue?',
            abort=True,
        )

    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance_name, NEO4J_PASSWORD)
        adapter.initialize_new_instance()
        click.echo(f"Database instance '{instance_name}' initialized successfully")
    except Exception as e:
        click.echo(f"Failed to initialize database instance: {str(e)}", err=True)


@db.command()
@click.argument("pdb_id", type=str)
@click.option(
    "--force", is_flag=True, help="Force upload even if structure already exists"
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload(
    pdb_id: str,
    force: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload a structure to the database using its PDB ID.

    PDB_ID should be a valid RCSB PDB identifier (e.g., '1j5e').
    """
    try:
        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Check if structure exists (unless force flag is used)
        if not force and adapter.check_structure_exists(pdb_id):
            click.echo(
                f"Structure {pdb_id.upper()} already exists in database. Use --force to overwrite."
            )
            return

        # Upload the structure
        click.echo(f"Uploading structure {pdb_id.upper()}...")
        adapter.add_total_structure(pdb_id, disable_exists_check=force)
        click.echo(f"Successfully uploaded structure {pdb_id.upper()}")

    except Exception as e:
        click.echo(f"Error uploading structure {pdb_id.upper()}: {str(e)}", err=True)
        raise click.Abort()


@db.command()
@click.option("--workers", "-w", type=int, default=4, help="Number of parallel workers")
@click.option(
    "--force", is_flag=True, help="Force upload even if structures already exist"
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload_all(
    workers: int,
    force: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload all available structures to the database in parallel.

    Uses GlobalOps to get list of structures and uploads them using multiple workers.
    """

    async def upload_structure(
        adapter: Neo4jAdapter, pdb_id: str, force: bool
    ) -> tuple[str, bool, Optional[str]]:
        """Upload a single structure and return result"""
        try:
            if not force and adapter.check_structure_exists(pdb_id):
                return pdb_id, False, "Already exists"

            adapter.add_total_structure(pdb_id, disable_exists_check=force)
            return pdb_id, True, None

        except Exception as e:
            return pdb_id, False, str(e)

    async def upload_worker(
        queue: asyncio.Queue, adapter: Neo4jAdapter, force: bool
    ) -> None:
        """Worker to process structures from the queue"""
        while True:
            try:
                pdb_id = await queue.get()
                result = await upload_structure(adapter, pdb_id, force)

                # Print result
                pdb_id, success, error = result
                if success:
                    click.echo(f"✓ {pdb_id}: Successfully uploaded")
                else:
                    click.echo(f"✗ {pdb_id}: {error}")

                queue.task_done()

            except asyncio.CancelledError:
                break

    async def main():
        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Get list of structures
        structures = GlobalOps.list_profiles()
        total = len(structures)

        click.echo(f"Found {total} structures to upload")

        # Create queue and workers
        queue = asyncio.Queue()
        worker_tasks = []

        # Start workers
        for _ in range(workers):
            task = asyncio.create_task(upload_worker(queue, adapter, force))
            worker_tasks.append(task)

        # Add structures to queue
        for pdb_id in structures:
            await queue.put(pdb_id)

        # Wait for all uploads to complete
        await queue.join()

        # Cancel workers
        for task in worker_tasks:
            task.cancel()
        await asyncio.gather(*worker_tasks, return_exceptions=True)

    try:
        # Run the async main function
        asyncio.run(main())
        click.echo("\nUpload complete!")

    except Exception as e:
        click.echo(f"Error during upload: {str(e)}", err=True)
        raise click.Abort()


@db.command()
@click.option("--workers", "-w", type=int, default=4, help="Number of parallel workers")
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show which structures would be uploaded without uploading",
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload_missing(
    workers: int,
    dry_run: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload structures that exist in RCSB but not in the local database.

    Uses GlobalOps.status_vs_rcsb() to identify missing structures and uploads them in parallel.
    """

    async def upload_structure(
        adapter: Neo4jAdapter, pdb_id: str
    ) -> tuple[str, bool, Optional[str]]:
        """Upload a single structure and return result"""
        try:
            adapter.add_total_structure(pdb_id, disable_exists_check=True)
            return pdb_id, True, None
        except Exception as e:
            return pdb_id, False, str(e)

    async def upload_worker(queue: asyncio.Queue, adapter: Neo4jAdapter) -> None:
        """Worker to process structures from the queue"""
        while True:
            try:
                pdb_id = await queue.get()
                result = await upload_structure(adapter, pdb_id)

                # Print result
                pdb_id, success, error = result
                if success:
                    click.echo(f"✓ {pdb_id}: Successfully uploaded")
                else:
                    click.echo(f"✗ {pdb_id}: {error}")

                queue.task_done()

            except asyncio.CancelledError:
                break

    async def main():
        # Get list of missing structures
        db_entries = Neo4jReader(
            Neo4jAdapter(
                uri=NEO4J_URI,
                user=NEO4J_USER,
                password=NEO4J_PASSWORD,
                current_db=NEO4J_CURRENTDB,
            )
        ).all_ids()
        missing_structures = GlobalOps.missing_db_entries(db_entries)
        total = len(missing_structures)

        if total == 0:
            click.echo("No missing structures found!")
            return

        for struct in missing_structures:
            click.echo(f"- {struct}")
        click.echo(
            f"Found {total} structures in RCSB that are not in the local database:"
        )

        if dry_run:
            click.echo("\nDry run - no structures were uploaded")
            return

        if not click.confirm(
            f"\nDo you want to proceed with uploading {total} structures?"
        ):
            click.echo("Upload cancelled")
            return

        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Create queue and workers
        queue = asyncio.Queue()
        worker_tasks = []

        # Start workers
        for _ in range(workers):
            task = asyncio.create_task(upload_worker(queue, adapter))
            worker_tasks.append(task)

        with click.progressbar(
            missing_structures, label="Queueing structures", length=total
        ) as structures:
            for pdb_id in structures:
                await queue.put(pdb_id)

        # Wait for all uploads to complete
        await queue.join()

        # Cancel workers
        for task in worker_tasks:
            task.cancel()
        await asyncio.gather(*worker_tasks, return_exceptions=True)

    try:
        # Run the async main function
        asyncio.run(main())
        click.echo("\nUpload complete!")

    except Exception as e:
        click.echo(f"Error during upload: {str(e)}", err=True)
        raise click.Abort()


@db.group()
@click.pass_context
def instance(ctx):
    """Manage Neo4j database instances"""
    pass


@instance.command()
@click.option(
    "--name", prompt="Database name", help="Name for the new database instance"
)
@click.option(
    "--initialize/--no-initialize",
    default=True,
    help="Initialize with constraints and basic data",
)
@click.option(
    "--force", "-f", is_flag=True, help="Force creation even if database exists"
)
def create(name: str, initialize: bool, force: bool):
    """Create a new Neo4j database instance.

    Examples:\n
    \b
    # Create a new database with prompts
    ribd.py db instance create

    \b
    # Create a specific database non-interactively
    ribd.py db instance create --name ribosome_test

    \b
    # Create without initialization
    ribd.py db instance create --name ribosome_test --no-initialize
    """
    try:
        # Connect to system database to manage instances
        adapter = Neo4jAdapter(
            NEO4J_URI,
            NEO4J_USER,
            "system",  # Use system database for management
            NEO4J_PASSWORD,
        )

        # Check if database exists
        if not force:
            with adapter.driver.session() as session:
                result = session.run("SHOW DATABASES WHERE name = $name", name=name)
                if result.single():
                    if not click.confirm(
                        f'Database "{name}" already exists. Override?'
                    ):
                        click.echo("Cancelled.")
                        return

        # Create database
        with adapter.driver.session() as session:
            session.run(f"CREATE DATABASE {name} IF NOT EXISTS")
            click.echo(f"Created database: {name}")

        if initialize:
            # Switch adapter to new database
            adapter = Neo4jAdapter(
                NEO4J_URI, NEO4J_USER, name, NEO4J_PASSWORD  # Use new database
            )

            with click.progressbar(
                length=3, label="Initializing database", show_pos=True
            ) as bar:
                # Initialize constraints
                adapter.init_constraints()
                bar.update(1)

                # Initialize polymer classes
                adapter.init_polymer_classes()
                bar.update(1)

                # Initialize phylogenies
                adapter.init_phylogenies()
                bar.update(1)

            click.echo("\nDatabase initialized successfully!")

    except Exception as e:
        click.echo(f"Error creating database: {str(e)}", err=True)
        raise


@instance.command()
def list():
    """List all Neo4j database instances"""
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, "system", NEO4J_PASSWORD)

        with adapter.driver.session() as session:
            result = session.run("SHOW DATABASES")
            databases = result.data()

            if not databases:
                click.echo("No databases found")
                return

            click.echo("\nAvailable databases:")
            for db in databases:
                status = "🟢" if db["currentStatus"] == "online" else "🔴"
                click.echo(f"{status} {db['name']:<20} ({db['currentStatus']})")

    except Exception as e:
        click.echo(f"Error listing databases: {str(e)}", err=True)


@instance.command(name="delete")
@click.argument("name", required=True)
@click.option("--force", "-f", is_flag=True, help="Force deletion without confirmation")
@click.pass_context
def delete_instance(ctx, name, force):
    """Delete a Neo4j database instance and all its data.

    Examples:\n
    \b
    # Delete an instance with confirmation
    ribd.py db instance delete ribosome2

    \b
    # Force delete without confirmation
    ribd.py db instance delete --force ribosome2
    """
    try:
        # Connect to system database for instance management
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, "system", NEO4J_PASSWORD)

        with adapter.driver.session() as session:
            # Check if instance exists
            result = session.run("SHOW DATABASES")
            databases = [record["name"] for record in result]

            if name not in databases:
                click.echo(f"Database instance '{name}' does not exist.", err=True)
                return

            # Check if it's the default database
            if name == "neo4j":
                click.echo("Cannot delete the default 'neo4j' database.", err=True)
                return

            # Get instance status
            result = session.run("SHOW DATABASE $name", name=name)
            status = result.single()["currentStatus"]

            if not force:
                click.echo(f"\nDatabase: {name}")
                click.echo(f"Status: {status}")
                click.confirm(
                    "Are you sure you want to delete this database and ALL its data?",
                    abort=True,
                )

            # Stop the database if it's running
            if status != "offline":
                click.echo("Stopping database...")
                session.run("STOP DATABASE $name", name=name)

            # Drop the database
            click.echo("Deleting database...")
            session.run("DROP DATABASE $name IF EXISTS", name=name)

            # Verify deletion
            result = session.run("SHOW DATABASES")
            databases = [record["name"] for record in result]

            if name not in databases:
                click.echo(f"Successfully deleted database instance '{name}'")
            else:
                click.echo(f"Failed to delete database instance '{name}'", err=True)

    except Exception as e:
        click.echo(f"Failed to delete database instance: {str(e)}", err=True)


if __name__ == "__main__":
    cli(obj={})
