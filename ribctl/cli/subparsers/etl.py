import asyncio
from pprint import pprint
import typing
import click
from click import Context

from ribctl.etl import etl_obtain
from ribctl.etl.etl_assets_ops import AssetClass, Assets

@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    if ctx.invoked_subcommand is None:

        stdin_text = click.get_text_stream('stdin')

        rcsb_id = ctx.obj['rcsb_id']
        # print(Assets.status_vs_rcsb())
        print(Assets(rcsb_id).status())

        # - whichassets are missing, 
        # - how up-to date with RCSB are we
        # - last error logs.
        # click.echo("Hi! You ran the 'etl' command without any subcommands.")

@etl.command()
@click.pass_context
@click.argument('assets', required=True,  nargs=-1,  type=click.Choice(AssetClass._member_names_))
@click.option('--overwrite', required=False,  type=bool)
def assets(ctx:Context, assets, overwrite):
    rcsb_id = ctx.obj['rcsb_id']
    print("Obtaining assets ", assets , " for {}.".format(rcsb_id)) 
    enums = list(map(AssetClass.from_str, assets))
    asyncio.run(etl_obtain.execute_asset_task_pool(etl_obtain.asset_routines(rcsb_id, enums, overwrite)))


