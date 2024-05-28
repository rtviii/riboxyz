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


    def list_ligands(self):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run("""//
                        match (n:Ligand)-[]-(r:RibosomeStructure) where not toLower(n.chemicalName) contains "ion" return properties(n), collect(r.rcsb_id)
                                        """).value()

            return session.execute_read(_)

    # TODO : {FETCH POLYMERS (RNA/PROTEIN)(BYSTRUCTURE/BYPOLYMERCLASS)}
    def list_polymers_filtered(self):
        ...

    def list_structs_filtered(self,
                            page: int,
                            search          : None| str                                           = None,
                            year            : None| typing.Tuple[int | None , int | None]         = None,
                            resolution      : None| typing.Tuple[float | None , float | None]     = None,
                            polymer_classes: None | list[PolynucleotideClass | PolypeptideClass ] = None,
                            source_taxa     : None| list[int]                                     = None,
                            host_taxa       : None| list[int]                                     = None ):

        query = """match (rib:RibosomeStructure)
with rib order by rib.rcsb_id desc\n""" + \
    \
( "\nwhere\n" if list(map(lambda x: x is not None, [search, year, resolution, polymer_classes, source_taxa, host_taxa]) ).count(True) > 0 else '' ) + \
\
( "toLower(rib.citation_title) + toLower(rib.rcsb_id) + toLower(rib.pdbx_keywords_text) + apoc.text.join(rib.citation_rcsb_authors, \"\")  contains '{}' \n".format(
    search) if search != '' else '' )  +\
( "{} ({} <= rib.citation_year and rib.citation_year <= {} or rib.citation_year is null)\n".format(
    "and" if search != '' else '', year[0] if year[0] is not None else 0, year[1] if year[1] is not None else 9999)  if year is not None else '') + \
( "{} {} < rib.resolution and rib.resolution < {}\n".format(
    "and" if search != '' or year !=None else '', resolution[0] if resolution[0] is not None else 0, resolution[1] if resolution[1] is not None else 9999) if resolution is not None else '') + \
( "{} ALL(x in [ {} ] where x in apoc.coll.flatten(collect{{ match (rib)-[]-(p:Polymer) return p.nomenclature }}))\n".format(
    "and" if search!= '' or year!=None or resolution != None else '',
    ', '.join(['"%s"' % w.value for w in polymer_classes])) if polymer_classes is not None else '' ) +\
( "{} ANY(tax in {} where tax in apoc.coll.flatten(collect{{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}}))\n".format(
    "and" if search!= '' or year!=None or resolution!=None or  polymer_classes!=None else '', source_taxa) if source_taxa is not None else '') + \
( "{} ANY(tax in {} where tax in apoc.coll.flatten(collect{{ match (rib)-[:host]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}}))\n".format(
    "and" if search!= '' or year!=None or resolution!=None or polymer_classes !=None or source_taxa!=None else '', host_taxa) if host_taxa is not None else '') + \
        \
"""
with collect(rib)[{}..{}] as rib, count(rib) as total_count 
unwind rib as ribosomes

optional match (l:Ligand)-[]-(ribosomes) 
with collect(PROPERTIES(l)) as ligands, ribosomes, total_count
match (rps:Protein)-[]-(ribosomes) 
with collect(PROPERTIES(rps)) as proteins, ligands, ribosomes, total_count
optional match (rna:RNA)-[]-(ribosomes) 
with collect(PROPERTIES(rna)) as rnas, proteins, ligands, ribosomes, total_count

with apoc.map.mergeList([{{proteins:proteins}},{{nonpolymeric_ligands:ligands}},{{rnas:rnas}},{{other_polymers:[]}}]) as rest, ribosomes, total_count
return collect(apoc.map.merge(ribosomes, rest)),  collect(distinct total_count)[0]
""".format((page-1)*20, page*20)

        print("=======Executing query:==========")
        print('\033[96m' + query + '\033[0m')
        print("===============================")

        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                return tx.run(query).values()
            return session.execute_read(_)

    # def list_structs_filtered(self,
    #                         search         : None|str=None,
    #                         year           : None| typing.Tuple[int | None , int | None] =None,
    #                         resolution     : None| typing.Tuple[float | None , float | None] =None, 
    #                         polymer_classes: None| list[PolynucleotideClass | PolypeptideClass ] =None,
    #                         source_taxa    : None| list[int] =None,
    #                         host_taxa      : None| list[int] =None ):
    #     with self.adapter.driver.session() as session:
    #         def _(tx: Transaction | ManagedTransaction):

    #             query = """//
    #                     match (rib:RibosomeStructure) 
    #                     with rib
    #                     order by rib.rcsb_id desc
    #                     where toLower(rib.citation_title) 
    #                           + toLower(rib.pdbx_keywords_text) 
    #                           + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
    #                     and rib.citation_year > 2020 
    #                     and rib.resolution < 3
    #                     and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
    #                     and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
    #                     with rib limit 10
    #                     optional match (l:Ligand)-[]-(rib) 
    #                     with collect(PROPERTIES(l)) as ligands, rib

    #                     match (rps:Protein)-[]-(rib) 
    #                     with collect(PROPERTIES(rps)) as proteins, ligands, rib

    #                     optional match (rna:RNA)-[]-(rib) 
    #                     with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

    #                     with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
    #                     return collect(apoc.map.merge(rib, rest))
    #                           """
    #             return tx.run(query).value()[0]
    #         return session.execute_read(_)

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






    # """
    # match (rib:RibosomeStructure) 
    # with rib
    # order by rib.rcsb_id desc 
    # where toLower(rib.citation_title) 
    #       + toLower(rib.pdbx_keywords_text) 
    #       + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
    # and rib.citation_year > 2020 
    # and rib.resolution < 3
    # and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
    # and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
    # with rib limit 10
    #     optional match (l:Ligand)-[]-(rib) 
    #     with collect(PROPERTIES(l)) as ligands, rib

    #     match (rps:Protein)-[]-(rib) 
    #     with collect(PROPERTIES(rps)) as proteins, ligands, rib

    #     optional match (rna:RNA)-[]-(rib) 
    #     with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

    # with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
    # return apoc.map.merge(rib, rest)
    # """
    # #* ................................I................
    # #* COUNT
    # """
    # match (rib:RibosomeStructure) 
    # with rib
    # order by rib.rcsb_id desc 
    # where toLower(rib.citation_title) 
    #       + toLower(rib.pdbx_keywords_text) 
    #       + apoc.text.join(rib.citation_rcsb_authors, "")  contains "complex" 
    # and rib.citation_year > 2020 
    # and rib.resolution < 3
    # and ALL(x in ["uL4", "uL22"] where x in apoc.coll.flatten(collect{match (rib)-[]-(p:Polymer) return p.nomenclature }) )
    # and ANY(tax in [9606] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}) )
    # return count(rib)
    # """

""""
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
                        
                        with collect(rib)[..10] as rib, count(rib) as total_count 
                        Unwind rib as ribs
                        optional match (l:Ligand)-[]-(ribs) 
                        with collect(PROPERTIES(l)) as ligands, ribs, total_count

                        match (rps:Protein)-[]-(ribs) 
                        with collect(PROPERTIES(rps)) as proteins, ligands, ribs, total_count

                        optional match (rna:RNA)-[]-(ribs) 
                        with collect(PROPERTIES(rna)) as rnas, proteins, ligands, ribs, total_count

                        with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, ribs, total_count
                        return collect(apoc.map.merge(ribs, rest)),  collect(distinct total_count)[0]
"""