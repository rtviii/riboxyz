
import typing
from neo4j import ManagedTransaction, Transaction
from neo4j_adapter.adapter import Neo4jAdapter

class Neo4jQuery():

    adapter: Neo4jAdapter 
    def __init__(self) -> None:
        self.adapter = Neo4jAdapter('bolt://localhost:7687', 'neo4j')
        pass

    def list_structs(self, filters=None, limit=None, offset=None):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//

        match (rib:RibosomeStructure) 
        with rib order by rib.rcsb_id limit 12

        optional match (l:Ligand)-[]-(rib) 
        with collect(PROPERTIES(l)) as ligands, rib

        match (rps:Protein)-[]-(rib) 
        with collect(PROPERTIES(rps)) as proteins, ligands, rib

        optional match (rna:RNA)-[]-(rib) 
        with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib
        
        with  apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
        return apoc.map.merge(rib, rest)
                                        """).value()

            return session.execute_read(_)

    def get_taxa(self, src_host:typing.Literal['source', 'host'])-> list[int]:
        def _(tx: Transaction | ManagedTransaction):
            return tx.run("""//
                    match (p:PhylogenyNode)-[k:{}]->(s:RibosomeStructure)
                    return collect(distinct(p.ncbi_tax_id))
            """.format(src_host)).single().value()

        with self.adapter.driver.session() as session:
            return session.execute_read(_)


dbqueries = Neo4jQuery()