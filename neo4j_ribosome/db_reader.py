from pprint import pprint
import typing
import sys

from ninja import Schema

from ribctl.lib.schema.types_ribosome import PolymerClass, PolynucleotideClass, PolypeptideClass
sys.dont_write_bytecode = True
from neo4j import ManagedTransaction, Transaction
from neo4j_ribosome.db_builder import Neo4jBuilder

"""This is the primary interface to the Neo4j instance. It queries the database and passes the results to the [ Django ] API.
DO NOT put validation logic/schema here. This is a pure interface to the database. All the conversions are done in the API layer.
"""

class FiltersSchema():
    search         : str
    year           : typing.Tuple[int | None , int | None]
    resolution     : typing.Tuple[float | None , float | None]
    polymer_classes: list[PolynucleotideClass | PolypeptideClass ]
    source_taxa    : list[int]
    host_taxa      : list[int]
    def __init__(self) -> None:
        pass


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

    def list_structs_filtered(self,
                            search         : None|str=None,
                            year           : None| typing.Tuple[int | None , int | None] =None,
                            resolution     : None| typing.Tuple[float | None , float | None] =None, 
                            polymer_classes: None| list[PolynucleotideClass | PolypeptideClass ] =None,
                            source_taxa    : None| list[int] =None,
                            host_taxa      : None| list[int] =None
 ):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):

                query = """//
                        match (rib:RibosomeStructure) 
                        with rib
                        order by rib.rcsb_id desc\n""" + ( """
                        where toLower(rib.citation_title) 
                              + toLower(rib.pdbx_keywords_text) 
                              + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
                        and rib.citation_year > 2020 
                        and rib.resolution < 3
                        and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
                        and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
                        """ if search is not None else '' ) + """
                        with rib limit 200
                        optional match (l:Ligand)-[]-(rib) 
                        with collect(PROPERTIES(l)) as ligands, rib

                        match (rps:Protein)-[]-(rib) 
                        with collect(PROPERTIES(rps)) as proteins, ligands, rib

                        optional match (rna:RNA)-[]-(rib) 
                        with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

                        with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
                        return collect(apoc.map.merge(rib, rest))
                              """
                pprint(query)
                return tx.run(query).value()[0]

            return session.execute_read(_)

    """
    match (rib:RibosomeStructure) 
    with rib
    order by rib.rcsb_id desc 
    where toLower(rib.citation_title) 
          + toLower(rib.pdbx_keywords_text) 
          + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
    and rib.citation_year > 2020 
    and rib.resolution < 3
    and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
    and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
    with rib limit 10
        optional match (l:Ligand)-[]-(rib) 
        with collect(PROPERTIES(l)) as ligands, rib

        match (rps:Protein)-[]-(rib) 
        with collect(PROPERTIES(rps)) as proteins, ligands, rib

        optional match (rna:RNA)-[]-(rib) 
        with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

    with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
    return apoc.map.merge(rib, rest)
    """
    #* ................................I................
    #* COUNT
    """
    match (rib:RibosomeStructure) 
    with rib
    order by rib.rcsb_id desc 
    where toLower(rib.citation_title) 
          + toLower(rib.pdbx_keywords_text) 
          + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
    and rib.citation_year > 2020 
    and rib.resolution < 3
    and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
    and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
    return count(rib)
    """


   




    # def list_struct_count(self, filters=None, limit=None, offset=None):
    #     with self.adapter.driver.session() as session:
    #         def _(tx: Transaction | ManagedTransaction):
    #             return tx.run("""//

    #     match (rib:RibosomeStructure) 
    #     with rib order by rib.rcsb_id desc limit 10

    #     optional match (l:Ligand)-[]-(rib) 
    #     with collect(PROPERTIES(l)) as ligands, rib

    #     match (rps:Protein)-[]-(rib) 
    #     with collect(PROPERTIES(rps)) as proteins, ligands, rib

    #     optional match (rna:RNA)-[]-(rib) 
    #     with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib
        
    #     with  apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
    #     return apoc.map.merge(rib, rest)
    #                                     """).value()

    #         return session.execute_read(_)



    def get_taxa(self, src_host:typing.Literal['source', 'host'])-> list[int]:
        def _(tx: Transaction | ManagedTransaction):
            return tx.run("""//
                    match (p:PhylogenyNode)-[k:{}]->(s:RibosomeStructure)
                    return collect(distinct(p.ncbi_tax_id))
            """.format(src_host)).single().value()

        with self.adapter.driver.session() as session:
            return session.execute_read(_)

dbqueries = Neo4jQuery()