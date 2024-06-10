import click
import sys
sys.dont_write_bytecode = True
from ribctl.cli.subparsers.ls import ls
from subparsers.etl import etl
from subparsers.db import db

@click.group()
@click.option('--debug/--no-debug', default=False, help='Enable/disable debug mode')
@click.option('--config', type=click.Path(exists=True), help='Path to the configuration file')
@click.option('--rcsb_id', required=False, help='RCSB ID')
@click.pass_context
def cli(ctx, debug, config, rcsb_id):
    if debug:
        click.echo('Debug mode is on')
    if config:
        click.echo(f'Using configuration file: {config}')
    ctx.ensure_object(dict)
    ctx.obj['rcsb_id'] = rcsb_id

cli.add_command(etl)
cli.add_command(db)
cli.add_command(ls)

@cli.command()
def init():
    """
    Initialize the application.
    """
    click.echo('Initializing the application...')

if __name__ == '__main__':
    cli(obj={})