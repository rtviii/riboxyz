import asyncio
from concurrent.futures import ALL_COMPLETED, Future, ThreadPoolExecutor, wait
from functools import partial
from pprint import pprint
import sys
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_reader import Neo4jReader
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from ribctl.global_ops import GlobalOps
from ribctl.ribosome_ops import StructureAssets, RibosomeOps, Structure

sys.dont_write_bytecode = True


# * Recipe for initializing a new instance from the RIBETL_DATA pool
# * - assumes the RibosomeStrucutre profiles are rendered
def full_upload():
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
    adapter.initialize_new_instance()
    futures: list[Future] = []

    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(StructureAssets.list_all_structs()):
            fut = executor.submit(partial(adapter.add_total_structure, rcsb_id, True))
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)


# * Add only the profiles that are missing versus the RCSB.
def rcsb_sync():
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(StructureAssets.status_vs_rcsb()):
            fut = executor.submit(partial(adapter.add_total_structure, rcsb_id, True))
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)


# * Alter only the core structure nodes given the profile is already rendered. (Does not alter the polymer nodes, etc.)
def upsert_all_structures():
    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
    futures: list[Future] = []
    with ThreadPoolExecutor(max_workers=10) as executor:
        for rcsb_id in sorted(GlobalView.list_all_structs()):
            fut = executor.submit(partial(adapter.upsert_structure_node, rcsb_id))
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)


def upsert_all_ligands():
    unique = {}
    for rcsb_id in sorted(StructureAssets.list_all_structs()):
        profile = RibosomeOps(rcsb_id).profile()
        for ligand in profile.nonpolymeric_ligands:
            if (not "ion" in ligand.chemicalName.lower()) and (
                ligand.chemicalId not in unique
            ):
                unique[ligand.chemicalId] = ligand

    adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)
    futures: list[Future] = []

    with ThreadPoolExecutor(max_workers=10) as executor:
        for ligand in list(unique.values()):
            fut = executor.submit(partial(adapter.upsert_ligand_node, ligand))
            futures.append(fut)
