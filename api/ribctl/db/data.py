from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import Protein
from api.schema.v0 import NeoStruct
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass


class QueryOps(Neo4jDB):
    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

    def get_all_structures(self)->list[NeoStruct]:

        with self.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return   tx.run("""//
                        match (ribs:RibosomeStructure) 
                                unwind ribs as struct

                                optional match (l:Ligand)-[]-(struct)
                                with collect(l.chemicalId) as ligands, struct

                                optional match (rps:Protein)-[]-(struct)
                                with ligands, struct, collect({
                                    auth_asym_id                   : rps.auth_asym_id,
                                    nomenclature                   : rps.nomenclature,
                                    entity_poly_seq_one_letter_code: rps.entity_poly_seq_one_letter_code
                                    }) as rps

                                optional match (struct_rnas:RNA)-[]-(struct)
                                with ligands, struct, rps, collect({
                                    auth_asym_id                   : struct_rnas.auth_asym_id,
                                    nomenclature                   : struct_rnas.nomenclature,
                                    entity_poly_seq_one_letter_code: struct_rnas.entity_poly_seq_one_letter_code
                                    }) as rnas

                                return struct, ligands,rps,rnas limit 3
                                        """).data()

            return session.read_transaction(_)

