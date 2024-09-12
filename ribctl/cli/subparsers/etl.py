import asyncio
import click
from click import Context
from loguru import logger
from ribctl.etl import etl_obtain
from ribctl.etl.etl_assets_ops import AssetClass, Assets
from ribctl.logs.loggers import get_etl_logger


@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    pass


@etl.command()
@click.pass_context
def all(ctx:Context):
        global_status = Assets.global_status()
        list(map(lambda a: a.name, list(AssetClass)))
        header = ( "[RCSB_ID \t" + "".join(["| {}".format(a.name) for a in list(AssetClass)]) + "]" )
        click.echo(header)
        for i, (k, v) in enumerate(global_status.items()):
            if i % 20 == 0:
                click.echo(header)
            assets_row = "{}\t\t".format(k) + "".join( ["| {}  ".format("X" if v[a] else " ") for a in list(AssetClass)] )
            click.echo(assets_row)


@etl.command()
@click.pass_context
@click.argument("assets", required=True, nargs=-1, type=click.Choice(AssetClass._member_names_) )
@click.option("--overwrite"  , required=False, is_flag=True, default=False)
@click.option("--rcsb_sync"  , required=False, is_flag=True, default=False)
@click.option("--all_structs", required=False, is_flag=True, default=False)
@click.option("--display", required=False, is_flag=True, default=False)
def assets(ctx: Context, assets, overwrite, rcsb_sync, all_structs, display):

    logger  = get_etl_logger()
    rcsb_id = ctx.obj['rcsb_id'].upper()
    assets  = list(map(AssetClass.from_str, assets))

    if rcsb_id is not None:

        routines = etl_obtain.asset_routines(rcsb_id, assets , overwrite)
        asyncio.run(etl_obtain.execute_asset_task_pool(routines))
        return 

    if rcsb_sync:
        for rcsb_id in Assets.status_vs_rcsb():
            print("RCSB Sync: Fetching assets for {}".format(rcsb_id))
            routines = etl_obtain.asset_routines(rcsb_id, [AssetClass.profile, AssetClass.cif] )
            asyncio.run(etl_obtain.execute_asset_task_pool(routines))
        return

    if all_structs:
        for rcsb_id in Assets.list_all_structs():
            print("[{}]".format(rcsb_id))
            try:
                routines = etl_obtain.asset_routines(rcsb_id, assets , overwrite)
                asyncio.run(etl_obtain.execute_asset_task_pool(routines))
                logger.info("Processed successfully {}: {}".format(rcsb_id, assets))
            except Exception as e:
                logger.error("Error processing {}: {}".format(rcsb_id, e))
                print("Error processing {}: {}".format(rcsb_id, e))
        return

    # asyncio.run(etl_obtain.execute_asset_task_pool())
