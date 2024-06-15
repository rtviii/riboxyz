from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
import os
import sys
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_reader import Neo4jQuery
from neo4j_ribosome.db_lib_builder import Neo4jBuilder
from ribctl.etl.etl_assets_ops import Assets, RibosomeOps, Structure

sys.dont_write_bytecode = True

# * Recipe for initializing a new instance from the RIBETL_DATA pool
# * - assumes the RibosomeStrucutre profiles are rendered
def full_upload():
    adapter = Neo4jBuilder(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
    adapter.initialize_new_instance()
    futures: list[Future] = []

    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(Assets.list_all_structs()):
            fut = executor.submit(partial(adapter.add_structure, rcsb_id, True))
            futures.append(fut)
        wait(futures, return_when=ALL_COMPLETED)

full_upload()