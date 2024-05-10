
from neo4j import ManagedTransaction, Transaction
from neo4j_adapter.adapter import Neo4jAdapter

class DBQuery():

    adapter: Neo4jAdapter 
    def __init__(self) -> None:
        self.adapter = Neo4jAdapter('bolt://localhost:7687', 'neo4j')
        pass

    def list_structs(self, filters=None, limit=None, offset=None):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//

        match (rib:RibosomeStructure) 
        with rib order by rib.rcsb_id limit 1

        optional match (l:Ligand)-[]-(rib) 
        with collect(PROPERTIES(l)) as ligands, rib

        match (rps:Protein)-[]-(rib) 
        with collect(PROPERTIES(rps)) as proteins, ligands, rib

        optional match (rna:RNA)-[]-(rib) 
        with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib
        
        with  apoc.map.mergeList([{proteins:proteins},{ligands:ligands},{rnas:rnas}]) as rest, rib
        return apoc.map.merge(rib, rest)
                                        """).data()

            return session.execute_read(_)


dbqueries = DBQuery()