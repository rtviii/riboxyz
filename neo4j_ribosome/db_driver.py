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

# def connect_all_structures_to_phylogenies():
#     adapter = Neo4jBuilder(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
#     for rib in Assets.list_all_structs():
#         try:
#             adapter.link_structure_to_phylogeny(rib)
#         except Exception as e:
#             print(e)


# print(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
# adapter = Neo4jBuilder(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)

# adapter.link_structure_to_phylogeny("8OVE")
# connect_all_structures_to_phylogenies()
# adapter = Neo4jBuilder('bolt://localhost:7687', 'neo4j')
# p = RibosomeAssets('8OVE').profile()
# for tax in [ *p.host_organism_ids, *p.src_organism_ids ]:
#     print(tax)
# td = 679895
# print(Taxid.get_lineage(td))
# print(Taxid.rank(td))
#     adapter._create_lineage(tax)
# adapter.link_structure_to_phylogeny('8OVE')
