from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
import sys
from neo4j_ribosome.db_reader import Neo4jQuery
from neo4j_ribosome.db_builder import Neo4jBuilder
from ribctl.etl.ribosome_assets import RibosomeAssets
sys.dont_write_bytecode = True
from dotenv import load_dotenv
load_dotenv('.env')



#* Recipe for initializing a new instance from the RIBETL_DATA pool
#* - assumes the RiboosomeStrucutre profiles are rendered

def full_upload(constrains:bool=True):
    #TODO :  ------- SANITY CHECKS
    adapter = Neo4jBuilder('bolt://localhost:7687', 'neo4j')
    if constrains:
        adapter.initialize_new_instance()

    futures:list[Future] =  []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in RibosomeAssets.list_all_structs():
            fut          = executor.submit(partial(adapter.add_structure, rcsb_id, False))
            futures.append(fut)
        wait(futures, return_when=ALL_COMPLETED)


# full_upload()


adapter = Neo4jBuilder('bolt://localhost:7687', 'neo4j')
for rib in RibosomeAssets.list_all_structs():
    print("Ligands of {}".format(rib))
    for l in RibosomeAssets(rib).profile().nonpolymeric_ligands:
        adapter.upsert_ligand_node(l, rib)
