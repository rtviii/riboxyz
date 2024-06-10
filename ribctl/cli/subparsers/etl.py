from pprint import pprint
import click
from click import Context

@click.group(invoke_without_command=True)
@click.pass_context
def etl(ctx: Context):
    if ctx.invoked_subcommand is None:
        click.echo("Hi! You ran the 'etl' command without any subcommands.")
    else:
        click.echo("hlelo with subcommand")


@etl.command()
@click.pass_context
def obtain(ctx):
    pprint(ctx)
    rcsb_id   = ctx.obj['rcsb_id']
    click.echo(f'Extracting {rcsb_id}')

@etl.command()
@click.pass_context
def transform(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Transforming {rcsb_id} ')
