import click

@click.group()
@click.option('--data_type', required=True, type=click.Choice(['text', 'binary', 'image']), help='Type of data')
@click.pass_context
def etl(ctx, data_type):
    """
    Extract, Transform, and Load (ETL) operations.
    """
    ctx.obj['data_type'] = data_type

@etl.command()
@click.pass_context
def extract(ctx):
    """
    Extract data from a source.
    """
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Extracting {rcsb_id} ({data_type})...')

@etl.command()
@click.pass_context
def transform(ctx):
    """
    Transform data.
    """
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Transforming {rcsb_id} ({data_type})...')

@etl.command()
@click.pass_context
def load(ctx):
    """
    Load data into a destination.
    """
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Loading {rcsb_id} ({data_type})...')