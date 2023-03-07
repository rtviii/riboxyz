from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import Protein
from api.schema.v0 import NeoStruct
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass


class QueryOps(Neo4jDB):

    def get_all_structures(self)->list[NeoStruct]:
        with self.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                match (n:RibosomeStructure) return n.rcsb_id
                """).values()
        structs =  session.read_transaction(_)

