from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
import os
import sys
from neo4j_ribosome import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_reader import Neo4jQuery
from neo4j_ribosome.db_lib_builder import Neo4jBuilder
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.libtax import Taxid

sys.dont_write_bytecode = True
from dotenv import load_dotenv

load_dotenv(".env")

# * Recipe for initializing a new instance from the RIBETL_DATA pool
# * - assumes the RiboosomeStrucutre profiles are rendered


def full_upload(constrains: bool = True):

    adapter = Neo4jBuilder("bolt://localhost:7687", "neo4j")
    if constrains:
        adapter.initialize_new_instance()

    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in RibosomeOps.list_all_structs():
            fut = executor.submit(partial(adapter.add_structure, rcsb_id, False))
            futures.append(fut)
        wait(futures, return_when=ALL_COMPLETED)

adapter = Neo4jBuilder(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)

def connect_all_structures_to_phylogenies():
    adapter = Neo4jBuilder("bolt://localhost:7687", "neo4j")
    for rib in RibosomeOps.list_all_structs():
        p = RibosomeOps(rib).profile()
        for tax in [*p.host_organism_ids, *p.src_organism_ids]:
            adapter._create_lineage(tax)
        adapter.link_structure_to_phylogeny(rib)

connect_all_structures_to_phylogenies()

# adapter = Neo4jBuilder('bolt://localhost:7687', 'neo4j')
# p = RibosomeAssets('8OVE').profile()
# for tax in [ *p.host_organism_ids, *p.src_organism_ids ]:
#     print(tax)
# td = 679895
# print(Taxid.get_lineage(td))
# print(Taxid.rank(td))
#     adapter._create_lineage(tax)
# adapter.link_structure_to_phylogeny('8OVE')
