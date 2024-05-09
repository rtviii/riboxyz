from neo4j import ManagedTransaction, Transaction
from ninja import NinjaAPI
from neo4j_adapter.adapter import Neo4jAdapter
from routers.router_struct import structure_router
from routers.router_classes import classification_router



class DBQuery():


    adapter: Neo4jAdapter 
    def __init__(self) -> None:
        self.adapter = Neo4jAdapter('bolt://localhost:7687', 'neo4j')
        pass
    def list_structs(self):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                        match (rib:RibosomeStructure) 
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

                                with ligands, rps, rnas, keys(struct) as keys, struct 
                                return struct, ligands,rps,rnas
                                        """).data()

            return session.execute_read(_)


dbqueries = DBQuery()
root_api = NinjaAPI(docs_url='/', openapi_url='/openapi', title='ribxz API', description='A programming interface for ribosome data.', version='0.1.0')
root_api.add_router('/structure', structure_router)
root_api.add_router('/polymer_class', classification_router)
