from pprint import pprint
import typing
import click
from click import Context

from ribctl.etl.ribosome_assets import Asset, AssetList 

@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    if ctx.invoked_subcommand is None:
        click.echo("Hi! You ran the 'etl' command without any subcommands.")
    else:
        click.echo("invoked with subcommand")

@etl.command()
@click.argument('data_type', required=True,  nargs=-1,  type=click.Choice(typing.get_args(Asset)))
@click.pass_context
def assets(ctx:Context, data_type):
    print(data_type)
    # debug = ctx.obj['DEBUG']
    # pprint(ctx)
    # click.echo(f'Debug mode: {debug}')
    # click.echo(f'Extracting data of type: {data_type}')
    # rcsb_id = ctx.obj['rcsb_id']
    # al      = AssetList([ 'profile', 'ligands' ])
    # print(al)
    # click.echo(f'Transforming {rcsb_id} ')

# @etl.command()
# @click.pass_context
# def obtain(ctx):
#     pprint(ctx)
#     rcsb_id   = ctx.obj['rcsb_id']
#     click.echo(f'Extracting {rcsb_id}')

# @etl.command()
# @click.pass_context
# def transform(ctx):
#     rcsb_id = ctx.obj['rcsb_id']
#     click.echo(f'Transforming {rcsb_id} ')

