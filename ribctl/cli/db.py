from functools import partial
import sys

sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.asset_manager.assets_structure import StructureAssets
from ribctl.global_ops import GlobalOps
from ribctl.ribosome_ops import RibosomeOps
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from ribctl import RIBETL_DATA
from ribctl.asset_manager.parallel_acquisition import process_chunk
import click
from typing import List, Tuple
import multiprocessing
from concurrent.futures import (
    ALL_COMPLETED,
    Future,
    ProcessPoolExecutor,
    ThreadPoolExecutor,
    wait,
)
from tqdm import tqdm
from ribctl.asset_manager.asset_types import AssetType


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
@click.argument("pdb_ids", type=PDBIDsParam(), nargs=-1)
@click.option("--workers", "-w", default=10, help="Number of worker threads")
@click.option(
    "--force", "-f", is_flag=True, help="Force upload even if structure exists"
)
@click.option(
    "--mode",
    "-m",
    type=click.Choice(["full", "structure", "ligands"], case_sensitive=False),
    default="full",
    help="Upload mode: full (all data), structure (only structure nodes), or ligands",
)
@click.option(
    "--instance", "-i", help="Neo4j database instance name", default=NEO4J_CURRENTDB
)
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
    all_pdb_ids = list(set(ctx.obj["piped_pdb_ids"] + sum(pdb_ids, [])))

    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return

    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance, NEO4J_PASSWORD)
    except Exception as e:
        click.echo(f"Failed to connect to database: {str(e)}", err=True)
        return

    click.echo(f"Starting upload to instance '{instance}'...")
    with click.progressbar(
        length=len(all_pdb_ids),
        label="Uploading structures",
        show_pos=True,
        file=sys.stderr,
    ) as bar:

        futures: list[Future] = []

        with ThreadPoolExecutor(max_workers=workers) as executor:
            for rcsb_id in sorted(all_pdb_ids):
                if mode == "full":
                    fut = executor.submit(
                        partial(adapter.add_total_structure, rcsb_id, force)
                    )
                elif mode == "structure":
                    fut = executor.submit(
                        partial(adapter.upsert_structure_node, rcsb_id)
                    )
                elif mode == "ligands":
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


@db.command()
@click.option(
    "--batch-size",
    "-b",
    default=100,
    help="Number of nodes to process before committing",
)
@click.option(
    "--sleep-time", "-s", default=2.0, help="Sleep time in seconds between batches"
)
@click.option(
    "--max-retries", "-r", default=3, help="Maximum number of retries per operation"
)
@click.option(
    "--force", "-f", is_flag=True, help="Force upload even if structure exists"
)
@click.option(
    "--mode",
    "-m",
    type=click.Choice(["full", "structure", "ligands"], case_sensitive=False),
    default="full",
    help="Upload mode: full (all data), structure (only structure nodes), or ligands",
)
@click.argument("pdb_ids", type=PDBIDsParam(), nargs=-1)
@click.pass_context
def upload(ctx, batch_size, sleep_time, max_retries, force, mode, pdb_ids):
    """Upload structures to Neo4j with memory-safe batching."""
    from time import sleep
    import random

    def process_with_retry(operation, retries):
        for attempt in range(retries):
            try:
                return operation()
            except Exception as e:
                if attempt == retries - 1:
                    raise e
                sleep_time = (attempt + 1) * 2 + random.random()
                sleep(sleep_time)

    all_pdb_ids = list(set(ctx.obj["piped_pdb_ids"] + sum(pdb_ids, [])))

    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return

    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)
    except Exception as e:
        click.echo(f"Failed to connect to database: {str(e)}", err=True)
        return

    # Split into batches
    batches = [
        all_pdb_ids[i : i + batch_size] for i in range(0, len(all_pdb_ids), batch_size)
    ]

    click.echo(f"Processing {len(all_pdb_ids)} structures in {len(batches)} batches")

    errors = []
    with tqdm(total=len(all_pdb_ids)) as pbar:
        for batch_idx, batch in enumerate(batches):
            try:
                # Process each structure in the batch
                for rcsb_id in batch:
                    try:

                        def upload_operation():
                            if mode == "full":
                                adapter.add_total_structure(rcsb_id, force)
                            elif mode == "structure":
                                adapter.upsert_structure_node(rcsb_id)
                            elif mode == "ligands":
                                profile = RibosomeOps(rcsb_id).profile
                                for ligand in profile.nonpolymeric_ligands:
                                    if not "ion" in ligand.chemicalName.lower():
                                        adapter.upsert_ligand_node(ligand, rcsb_id)

                        process_with_retry(upload_operation, max_retries)
                        pbar.update(1)
                        pbar.set_description(f"Processed {rcsb_id}")

                    except Exception as e:
                        error_msg = f"Error processing {rcsb_id}: {str(e)}"
                        errors.append((rcsb_id, error_msg))
                        click.echo(f"\n{error_msg}", err=True)
                        pbar.update(1)

                # Sleep between batches to allow garbage collection
                if batch_idx < len(batches) - 1:
                    sleep(sleep_time)

            except Exception as e:
                error_msg = f"Failed to process batch: {str(e)}"
                click.echo(f"\n{error_msg}", err=True)
                errors.extend((str(pdb_id), error_msg) for pdb_id in batch)

    # Final report
    total_processed = len(all_pdb_ids)
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


@etl.command(name="get_all")
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
def get_all(
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
    """Download or generate assets for all available structures.

    Requires explicit specification of which asset types to acquire using --asset-types/-t option.

    Examples:\n
    \b
    # Get specific assets for all structures
    ribd.py etl get_all --asset-types MMCIF STRUCTURE_PROFILE

    \b
    # Get assets with custom settings
    ribd.py etl get_all --asset-types MMCIF PTC --workers 4 --chunk-size 10

    \b
    # Process with retries and delays
    ribd.py etl get_all --asset-types STRUCTURE_PROFILE --delay 5 --max-retries 5
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
        all_structures = GlobalOps.list_profiles()
        if not all_structures:
            click.echo("No structures found to process", err=True)
            return

        click.echo(f"Found {len(all_structures)} structures to process")
    except Exception as e:
        click.echo(f"Failed to get structure list: {str(e)}", err=True)
        return

    # Convert asset type names to enums
    selected_asset_types = [AssetType[name] for name in asset_types]
    assets_str = ", ".join(ast.name for ast in selected_asset_types)

    # Inform user about the operation
    click.echo(f"Getting assets [{assets_str}] for {len(all_structures)} structures...")
    click.echo(f"Using {workers} workers, {chunk_size} structures per chunk")
    click.echo(f"Delay between chunks: {delay}s")

    # Split all structures into chunks
    chunks = [
        all_structures[i : i + chunk_size]
        for i in range(0, len(all_structures), chunk_size)
    ]

    errors = []
    with tqdm(total=len(all_structures)) as pbar:
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
    total_processed = len(all_structures)
    failed = len(errors)
    succeeded = total_processed - failed

    click.echo("\nProcessing Summary:")
    click.echo(f"Total structures processed: {total_processed}")
    click.echo(f"Successfully processed: {succeeded}")
    click.echo(f"Failed: {failed}")

    if errors:
        click.echo("\nErrors occurred during processing:")
        for rcsb_id, error in errors:
            click.echo(f"  - {rcsb_id}: {error}")
    else:
        click.echo("\nAll structures processed successfully")
    # Final report
    total_processed = len(all_structures)
    failed = len(errors)
    succeeded = total_processed - failed

    click.echo("\nProcessing Summary:")
    click.echo(f"Total structures processed: {total_processed}")
    click.echo(f"Successfully processed: {succeeded}")
    click.echo(f"Failed: {failed}")

    if errors:
        click.echo("\nErrors occurred during processing:")
        for rcsb_id, error in errors:
            click.echo(f"  - {rcsb_id}: {error}")
    else:
        click.echo("\nAll structures processed successfully")


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


@db.command(name="upload_all")
@click.option(
    "--force", "-f", is_flag=True, help="Force upload even if structures exist"
)
@click.option(
    "--mode",
    "-m",
    type=click.Choice(["full", "structure", "ligands"], case_sensitive=False),
    default="full",
    help="Upload mode: full (all data), structure (only structure nodes), or ligands",
)
@click.pass_context
def upload_all(ctx, force, mode):
    """Upload all available structures one at a time."""
    try:
        all_structures = sorted(GlobalOps.list_profiles())
        if not all_structures:
            click.echo("No structures found to upload", err=True)
            return

        click.echo(f"Found {len(all_structures)} structures to process")

        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD)

        # Process one at a time
        errors = []
        with tqdm(total=len(all_structures), desc="Uploading structures") as pbar:
            for rcsb_id in all_structures:
                try:
                    if mode == "full":
                        adapter.add_total_structure(rcsb_id, force)
                    elif mode == "structure":
                        adapter.upsert_structure_node(rcsb_id)
                    elif mode == "ligands":
                        profile = RibosomeOps(rcsb_id).profile
                        for ligand in profile.nonpolymeric_ligands:
                            if not "ion" in ligand.chemicalName.lower():
                                adapter.upsert_ligand_node(ligand, rcsb_id)

                    pbar.update(1)
                    pbar.set_description(f"Processed {rcsb_id}")

                except Exception as e:
                    error_msg = f"Error processing {rcsb_id}: {str(e)}"
                    errors.append((rcsb_id, error_msg))
                    click.echo(f"\nError with {rcsb_id}: {str(e)}", err=True)
                    pbar.update(1)
                    continue  # Move to next structure

        # Final report
        total = len(all_structures)
        failed = len(errors)
        succeeded = total - failed

        click.echo("\nUpload Summary:")
        click.echo(f"Total structures processed: {total}")
        click.echo(f"Successfully uploaded: {succeeded}")
        click.echo(f"Failed: {failed}")

        if errors:
            click.echo("\nErrors occurred during upload:")
            for rcsb_id, error in errors:
                click.echo(f"  - {rcsb_id}: {error}")
        else:
            click.echo("\nAll structures uploaded successfully")

    except Exception as e:
        click.echo(f"An error occurred: {str(e)}", err=True)
        raise
