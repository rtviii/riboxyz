import click

@click.group()
@click.option('--debug/--no-debug', default=False, help='Enable/disable debug mode')
@click.option('--config', type=click.Path(exists=True), help='Path to the configuration file')
@click.option('--rcsb_id', required=False, help='RCSB ID')

@click.pass_context
def cli(ctx, debug, config, rcsb_id):
    """
    This is the main CLI application.
    """
    if debug:
        click.echo('Debug mode is on')
    if config:
        click.echo(f'Using configuration file: {config}')

    ctx.ensure_object(dict)
    ctx.obj['rcsb_id'] = rcsb_id

@cli.command()
def init():
    """
    Initialize the application.
    """
    click.echo('Initializing the application...')

@cli.group()
@click.pass_context
@click.option('--data_type', required=True, type=click.Choice(['text', 'binary', 'image']), help='Type of data')
def etl(ctx):
    """
    Extract, Transform, and Load (ETL) operations.
    """
    rcsb_id = ctx.obj['rcsb_id']

@etl.command()
@click.pass_context
def extract(ctx):
    """
    Extract data from a source.
    """
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Extracting {rcsb_id}...')

@etl.command()
@click.pass_context
def transform(ctx):
    """
    Transform data.
    """
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Transforming {rcsb_id}...')

@etl.command()
@click.pass_context
def load(ctx):
    """
    Load data into a destination.
    """
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Loading {rcsb_id}...')

@cli.group()
@click.pass_context
def lig(ctx):
    """
    Linguistic operations.
    """
    rcsb_id = ctx.obj['rcsb_id']

@lig.command()
@click.pass_context
def tokenize(ctx):
    """
    Tokenize text.
    """
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Tokenizing {rcsb_id}...')

@lig.command()
@click.pass_context
def lemmatize(ctx):
    """
    Lemmatize text.
    """
    rcsb_id = ctx.obj['rcsb_id'li()