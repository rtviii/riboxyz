import click

@click.command()
@click.pass_context
def ls(ctx):
    rcsb_id = ctx.obj['rcsb_id']
    click.echo(f'Listing {rcsb_id}...')