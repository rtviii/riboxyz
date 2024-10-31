CHEMIDS = [
    "DI0",
    "PAR",
    "SRY",
    "DOL",
    "CLM",
    "ANM",
    "FUA",
    "SPS",
    "ZIT",
    "84G",
    "ACY",
    "TAC",
    "VIR",
    "CLY",
    "ERY",
    "AM2",
    "3QB",
    "TEL",
    "CAI",
    "AKN",
    "SCM",
    "62B",
    "GET",
    "TYK",
    "U7V",
    "PUY",
    "TAO",
    "MUL",
    "84D",
    "CTY",
    "V7A",
    "KAN",
    "EDS",
    "TOY",
    "EM1",
    "ZLD",
    "MAU",
    "G34",
]

import json
from pprint import pprint
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from neo4j_ribosome.db_lib_reader import dbqueries
from ribctl.etl.ribosome_ops import Assets
from ribctl.lib.libbsite import bsite_ligand


def all_structs_for_lig(chemid: str):
    adapter = Neo4jAdapter(
        uri="bolt://localhost:7687", user="neo4j", current_db="ribxz", password=""
    )
    Q = """match (n:Ligand {chemicalId:$chemid})-[]-(s:RibosomeStructure) return s.rcsb_id, s.src_organism_names, s.src_organism_ids"""
    with adapter.driver.session() as session:

        def _(tx):
            return tx.run(Q, {"chemid": chemid}).values()

        return session.execute_read(_)


#! save
for CHEMID in CHEMIDS:
    for i in all_structs_for_lig(CHEMID):
        rcsb_id = i[0]
        path    = Assets(rcsb_id).paths.binding_site(CHEMID)
        bsite   = bsite_ligand(CHEMID, rcsb_id)

        with open(path, "w+") as f:
            json.dump(bsite, f, indent=4)
            print("Saved: ", path)

        path2 = "/home/rtviii/dev/riboxyz/antibiotic_bsites/{}_{}.json".format( CHEMID, rcsb_id )
        with open(path2, "w+") as f:
            json.dump(bsite, f, indent=4)
            print("Saved: ", path2)

# subprocess.run(["chimerax", "--script", "/home/rtviii/dev/riboxyz/chimerax/ligand_vis.py", "--", CHEMID])

# x chimerax vis
