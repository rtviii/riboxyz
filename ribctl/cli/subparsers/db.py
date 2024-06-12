import click

@click.group()
@click.pass_context
def db(ctx):
    pass

@db.command()
@click.pass_context
def migrate(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Migrating {rcsb_id}...')

@db.command()
@click.pass_context
def seed(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Seeding {rcsb_id}...')