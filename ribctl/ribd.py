from pprint import pprint
import click
import sys
sys.dont_write_bytecode = True
sys.path.append('/home/rtviii/dev/riboxyz')
from ribctl.cli.subparsers.lig import lig
from cli.subparsers.ls import ls
from ribctl.cli.subparsers.etl import  etl
from cli.subparsers.db import db

ce = click.echo

@click.group()
@click.option('--debug/--no-debug', default=False, help='Enable/disable debug mode')
@click.option('--config', type=click.Path(exists=True), help='Path to the configuration file')
@click.option('--rcsb_id','-s', required=False, help='RCSB ID')
@click.pass_context
def ribd(ctx, debug, config, rcsb_id:str):
    if debug:
        click.echo('Debug mode is on')
    if config:
        click.echo(f'Using configuration file: {config}')

    ctx.ensure_object(dict)
    if rcsb_id is not None:
        ctx.obj['rcsb_id'] = rcsb_id.upper()

ribd.add_command(etl)
ribd.add_command(lig)
ribd.add_command(db)
ribd.add_command(ls)

@ribd.command()
def init():
    click.echo('Initializing the application...')

if __name__ == '__main__':
    ribd(obj={})