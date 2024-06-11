from pprint import pprint
import typing
import click
from click import Context

from ribctl.etl import obtain
from ribctl.etl.ribosome_assets import Asset

@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    if ctx.invoked_subcommand is None:
        click.echo("Hi! You ran the 'etl' command without any subcommands.")
    else:
        click.echo("invoked with subcommand")

@etl.command()
@click.pass_context
@click.argument('assets', required=True,  nargs=-1,  type=click.Choice(Asset._member_names_))
def assets(ctx:Context, assets):
    rcsb_id = ctx.obj['rcsb_id']

    print("Obtaining assets ", assets , " for {}.".format(rcsb_id)) 
    enums = list(map(Asset.from_str, assets))
    print("got assets",  enums)
    [print(x) for x in obtain.obtain_assets(rcsb_id, enums)]
# 

