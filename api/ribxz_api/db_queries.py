import typing
from neo4j import ManagedTransaction, Transaction
from neo4j_ribosome.adapter import Neo4jBuilder

"""This is the primary interface to the Neo4j instance. It queries the database and passes the results to the [ Django ] API.
DO NOT put validation logic/schema here. This is a pure interface to the database. All the conversions are done in the API layer.
"""

class Neo4jQuery():

    adapter: Neo4jBuilder 
    def __init__(self) -> None:
        self.adapter = Neo4jBuilder('bolt://localhost:7687', 'neo4j')
        pass

    def list_chains_by_struct(self, filters=None, limit=None, offset=None):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                # TODO : Fix the limit here (didn't have all structures taged with polymers in dev)
                return tx.run("""//
                        match (rib:RibosomeStructure) 
                        with rib order by rib.rcsb_id desc  limit 120000000 
                        match (poly:Polymer)-[]-(rib)
                        with collect({
                            nomenclature: poly.nomenclature,
                            auth_asym_id: poly.auth_asym_id,
                            entity_poly_polymer_type: poly.entity_poly_polymer_type,
                            entity_poly_seq_length:poly.entity_poly_seq_length
                        }) as polymers, rib
                        return {rcsb_id: rib.rcsb_id,polymers: polymers}
                                        """).value()

            return session.execute_read(_)

    def list_structs(self, filters=None, limit=None, offset=None):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//

        match (rib:RibosomeStructure) 
        with rib order by rib.rcsb_id limit 40

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