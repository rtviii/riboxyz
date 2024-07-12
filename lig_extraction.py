CHEMIDS = [ "DI0", "PAR", "SRY", "DOL", "CLM", "ANM", "FUA", "SPS",
"ZIT", "84G", "ACY", "TAC", "VIR", "CLY", "ERY", "AM2", "3QB", "TEL", "CAI", "AKN",
"SCM", "62B", "GET", "TYK", "U7V", "PUY", "TAO", "MUL", "84D", "CTY", "V7A", "KAN", "EDS", "TOY", "EM1", "ZLD", "MAU", "G34"]

from pprint import pprint
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from neo4j_ribosome.db_lib_reader import dbqueries

def all_structs_for_lig(chemid:str):
    adapter = Neo4jAdapter(uri="bolt://localhost:7687", user="neo4j", current_db="ribxz", password="")
    Q = """match (n:Ligand {chemicalId:"TEL"})-[]-(s:RibosomeStructure) return s.rcsb_id, s.src_organism_names, s.src_organism_ids"""
    with adapter.driver.session() as session:
            def _(tx):
                return tx.run(Q).values()
            return session.execute_read(_)

pprint(all_structs_for_lig("TEL"))
# extract nbhd
