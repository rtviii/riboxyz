import click

@click.group()
@click.option('--data_type', required=True, type=click.Choice(['text', 'binary', 'image']), help='Type of data')
@click.pass_context
def etl(ctx, data_type):
    ctx.obj['data_type'] = data_type

@etl.command()
@click.pass_context
def extract(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Extracting {rcsb_id} ({data_type})...')

@etl.command()
@click.pass_context
def transform(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Transforming {rcsb_id} ({data_type})...')

@etl.command()
@click.pass_context
def load(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    data_type = ctx.obj['data_type']
    click.echo(f'Loading {rcsb_id} ({data_type})...')