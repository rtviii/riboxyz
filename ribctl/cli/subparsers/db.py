import click

from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_driver import upsert_all_structures
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from neo4j_ribosome.db_lib_reader import Neo4jReader


adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
reader  = Neo4jReader(adapter)

@click.group()
@click.pass_context
def db(ctx):
    pass

@db.command()
@click.pass_context
@click.option("--all_structs", required=False, is_flag=True, default=False)
def upsert(ctx, all_structs):
    if all_structs:
        upsert_all_structures()
        return

    rcsb_id = ctx.obj['rcsb_id']
    if rcsb_id !=None:
        adapter.upsert_structure_node(rcsb_id)
        click.echo(f'Upserted {rcsb_id} from profile...')

@db.command()
@click.pass_context
@click.option("--node_types", required=False, is_flag=True, default=True)
def ls(ctx, node_types):
    if node_types:
        click.echo(reader.node_types())