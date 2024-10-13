from pprint import pprint
import typing
import sys
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from ribctl.lib.schema.types_ribosome import (
    PolynucleotideClass,
    PolynucleotideClass,
    PolypeptideClass,
)

sys.dont_write_bytecode = True
from neo4j import ManagedTransaction, Transaction
from neo4j_ribosome.db_lib_builder import Neo4jAdapter

"""This is the primary interface to the Neo4j instance. It queries the database and passes the results to the [ Django ] API.
DO NOT put validation logic/schema here. This is a pure interface to the database. All the conversions are done in the API layer.
"""

from pydantic import BaseModel, Field
from typing import Optional, List, Literal


class FilterParams(BaseModel):

    cursor: Optional[str] = None
    limit: int = Field(default=20, ge=1, le=100)
    year: Optional[tuple[Optional[int], Optional[int]]] = None
    search: Optional[str] = None
    resolution: Optional[tuple[Optional[float], Optional[float]]] = None
    polymer_classes: Optional[List[PolynucleotideClass | PolypeptideClass]] = None
    source_taxa: Optional[List[int]] = None
    host_taxa: Optional[List[int]] = None
    subunit_presence: Optional[Literal["SSU+LSU", "LSU", "SSU"]] = None


class Neo4jReader:

    adapter: Neo4jAdapter

    def __init__(self, adapter: Neo4jAdapter | None = None) -> None:
        if adapter:
            self.adapter = adapter
            return
        self.adapter = Neo4jAdapter(
            NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB, NEO4J_PASSWORD
        )

    def node_types(self):
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """match (n)return collect(distinct labels(n))"""
                ).value()[0]

            return session.execute_read(_)

    def tax_dict(self):
        """All taxonomic ids present in the database mapped to their scientific name per NCBI"""
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """match (t:PhylogenyNode) return collect([t.ncbi_tax_id,t.scientific_name])"""
                ).value()[0]

            return session.execute_read(_)

    def all_ids(self):
        """All taxonomic ids present in the database mapped to their scientific name per NCBI"""
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """match (r:RibosomeStructure) return collect(r.rcsb_id)"""
                ).value()[0]

            return session.execute_read(_)

    def polymer_classes_stats(self):
        """All taxonomic ids present in the database mapped to their scientific name per NCBI"""
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """match (n:PolymerClass)-[]-(p:Polymer) with  count(p) as c, n return collect([n.class_id, c])"""
                ).value()[0]

            return session.execute_read(_)

    def random_structure(self):
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """
match (rib: RibosomeStructure) 
with rand() as r, rib
order by r limit 1

optional match (l:Ligand)-[]-(rib) 
with collect(PROPERTIES(l)) as ligands, rib
match (rps:Protein)-[]-(rib) 
with collect(PROPERTIES(rps)) as proteins, ligands, rib
optional match (rna:RNA)-[]-(rib) 
with collect(PROPERTIES(rna)) as rnas, proteins, ligands, rib

with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, rib
return apoc.map.merge(rib, rest)
        """
                ).value()

            return session.execute_read(_)

    def list_chains_by_struct(self, filters=None, limit=None, offset=None):
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """//
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
                                        """
                ).value()

            return session.execute_read(_)

    def structures_overview(self):
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(
                    """//
    match (n:RibosomeStructure)-[:source]-(p:PhylogenyNode) 
    return collect({rcsb_id:n.rcsb_id, tax_id: p.ncbi_tax_id, tax_name:p.scientific_name, mitochondrial:n.mitochondrial, title:n.citation_title})"""
                ).value()[0]

            return session.execute_read(_)

    def list_ligands(self, nodes_only: bool = False):
        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                if nodes_only:
                    return tx.run(
                        """match (l:Ligand) where not toLower(l.chemicalName) contains "ion" with l return collect(properties(l)) """
                    ).values()
                else:
                    return tx.run(
                        """match (l:Ligand)-[]-(r:RibosomeStructure) where not toLower(l.chemicalName) contains "ion" 
                                  with l, r 
                                  match (r)-[:source]-(p:PhylogenyNode) 
                                  return properties(l), collect({rcsb_id: r.rcsb_id , tax_node: properties(p)})"""
                    ).values()

            return session.execute_read(_)

    def list_polymers_filtered_by_polymer_class(
        self, page: int, polymer_class: PolynucleotideClass
    ):

        query_by_polymer_class = """
         match (poly:Polymer)-[]-(pc:PolymerClass {{class_id:"{}" }})
         with poly order by poly.nomenclature desc
         with collect(poly)[{}..{}] as poly, count(poly) as pcount
         unwind poly as polys
         return collect(properties(polys)), pcount
        """.format(
            polymer_class.value, (page - 1) * 50, page * 50
        )

        print("=======Executing polymer query by polymer class:==========")
        print("\033[96m" + query_by_polymer_class + "\033[0m")
        print("===============================")

        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(query_by_polymer_class).values()

            return session.execute_read(_)

    def list_polymers_filtered_by_structure(
        self,
        page: int,
        search: None | str = None,
        year: None | typing.Tuple[int | None, int | None] = None,
        resolution: None | typing.Tuple[float | None, float | None] = None,
        polymer_classes: None | list[PolynucleotideClass | PolypeptideClass] = None,
        source_taxa: None | list[int] = None,
        host_taxa: None | list[int] = None,
    ):

        query_by_structure = (
            """match (rib:RibosomeStructure)
with rib order by rib.rcsb_id desc\n"""
            + (
                "\nwhere\n"
                if list(
                    map(
                        lambda x: x is not None,
                        [
                            search,
                            year,
                            resolution,
                            polymer_classes,
                            source_taxa,
                            host_taxa,
                        ],
                    )
                ).count(True)
                > 0
                else ""
            )
            + (
                "toLower(rib.citation_title) + toLower(rib.rcsb_id) + toLower(rib.pdbx_keywords_text) + apoc.text.join(rib.citation_rcsb_authors, \"\")  contains '{}' \n".format(
                    search
                )
                if search != ""
                else ""
            )
            + (
                "{} ({} <= rib.citation_year and rib.citation_year <= {} or rib.citation_year is null)\n".format(
                    "and" if search != "" else "",
                    year[0] if year[0] is not None else 0,
                    year[1] if year[1] is not None else 9999,
                )
                if year is not None
                else ""
            )
            + (
                "{} {} < rib.resolution and rib.resolution < {}\n".format(
                    "and" if search != "" or year != None else "",
                    resolution[0] if resolution[0] is not None else 0,
                    resolution[1] if resolution[1] is not None else 9999,
                )
                if resolution is not None
                else ""
            )
            + (
                "{} ALL(x in [ {} ] where x in apoc.coll.flatten(collect{{ match (rib)-[]-(p:Polymer) return p.nomenclature }}))\n".format(
                    "and" if search != "" or year != None or resolution != None else "",
                    ", ".join(['"%s"' % w.value for w in polymer_classes]),
                )
                if polymer_classes is not None
                else ""
            )
            # # !------------- TAXONOMY
            #             + (
            #                 "{} ANY(tax in {} where tax in apoc.coll.flatten(collect{{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}}))\n".format(
            #                     (
            #                         "and"
            #                         if search != ""
            #                         or year != None
            #                         or resolution != None
            #                         or polymer_classes != None
            #                         else ""
            #                     ),
            #                     source_taxa,
            #                 )
            #                 if source_taxa is not None
            #                 else ""
            #             )
            #             + (
            #                 "{} ANY(tax in {} where tax in apoc.coll.flatten(collect{{ match (rib)-[:host]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}}))\n".format(
            #                     (
            #                         "and"
            #                         if search != ""
            #                         or year != None
            #                         or resolution != None
            #                         or polymer_classes != None
            #                         or source_taxa != None
            #                         else ""
            #                     ),
            #                     host_taxa,
            #                 )
            #                 if host_taxa is not None
            #                 else ""
            #             )
            # !------------- TAXONOMY NEW
            + (
                "{} exists{{ MATCH (rib)-[:belongs_to_lineage_source]-(p:PhylogenyNode ) where p.ncbi_tax_id in {} }}\n".format(
                    (
                        "and"
                        if search != ""
                        or year != None
                        or resolution != None
                        or polymer_classes != None
                        else ""
                    ),
                    source_taxa,
                )
                if source_taxa is not None
                else ""
            )
            + (
                "{} exists{{ MATCH (rib)-[:belongs_to_lineage_host]-(p:PhylogenyNode ) where p.ncbi_tax_id in {} }}\n".format(
                    (
                        "and"
                        if search != ""
                        or year != None
                        or resolution != None
                        or polymer_classes != None
                        else ""
                    ),
                    host_taxa,
                )
                if host_taxa is not None
                else ""
            )
            # !------------- TAXONOMY NEW
            + """
 with collect(rib) as rib
 unwind rib as ribosomes

 match (poly:Polymer)-[]-(ribosomes)
 with poly order by poly.nomenclature desc
 with collect(poly)[{}..{}] as poly, count(poly) as pcount
 unwind poly as polys
 return collect(properties(polys)), pcount
""".format(
                (page - 1) * 50, page * 50
            )
        )

        print("=======Executing polumer query by structure:==========")
        print("\033[96m" + query_by_structure + "\033[0m")
        print("===============================")

        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                return tx.run(query_by_structure).values()

            return session.execute_read(_)

    def list_structs_filtered(self, filters: FilterParams):
        query_parts = ["MATCH (rib:RibosomeStructure)"]
        where_clauses = []
        params = {"limit": filters.limit}

        if filters.cursor:
            where_clauses.append("rib.rcsb_id < $cursor")
            params["cursor"] = filters.cursor

        if filters.search:
            where_clauses.append(
                """
                toLower(rib.citation_title) + toLower(rib.rcsb_id) + 
                toLower(rib.pdbx_keywords_text) + 
                toLower(reduce(acc = '', str IN rib.src_organism_names | acc + str)) + 
                toLower(reduce(acc = '', str IN rib.host_organism_names | acc + str)) +
                toLower(reduce(acc = '', str IN rib.citation_rcsb_authors | acc + str))
                CONTAINS $search
            """
            )
            params["search"] = filters.search.lower()

        if filters.year:
            start, end = filters.year
            if start is not None:
                where_clauses.append(
                    "(datetime(rib.deposition_date).year >= $year_start AND rib.deposition_date IS NOT NULL)"
                )
                params["year_start"] = start
            if end is not None:
                where_clauses.append(
                    "(datetime(rib.deposition_date).year <= $year_end AND rib.deposition_date IS NOT NULL)"
                )
                params["year_end"] = end

        if filters.resolution:
            start, end = filters.resolution
            if start is not None:
                where_clauses.append("rib.resolution > $resolution_start")
                params["resolution_start"] = float(start)
            if end is not None:
                where_clauses.append("rib.resolution < $resolution_end")
                params["resolution_end"] = float(end)

        if filters.polymer_classes:
            where_clauses.append(
                "ALL(x IN $polymer_classes WHERE x IN apoc.coll.flatten(collect{ MATCH (rib)-[]-(p:Polymer) RETURN p.nomenclature }))"
            )
            params["polymer_classes"] = [pc.value for pc in filters.polymer_classes]

        if filters.source_taxa:
            where_clauses.append(
                "EXISTS{ MATCH (rib)-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }"
            )
            params["source_taxa"] = filters.source_taxa

        if filters.host_taxa:
            where_clauses.append(
                "EXISTS{ MATCH (rib)-[:belongs_to_lineage_host]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $host_taxa }"
            )
            params["host_taxa"] = filters.host_taxa

        if filters.subunit_presence:
            if filters.subunit_presence == "SSU+LSU":
                where_clauses.append(
                    '"lsu" IN rib.subunit_presence AND "ssu" IN rib.subunit_presence'
                )
            elif filters.subunit_presence == "LSU":
                where_clauses.append(
                    '"lsu" IN rib.subunit_presence AND NOT "ssu" IN rib.subunit_presence'
                )
            elif filters.subunit_presence == "SSU":
                where_clauses.append(
                    '"ssu" IN rib.subunit_presence AND NOT "lsu" IN rib.subunit_presence'
                )

        if where_clauses:
            query_parts.append("WHERE " + " AND ".join(where_clauses))

        query_parts.extend(
            [
                "WITH rib",
                "ORDER BY rib.rcsb_id DESC",
                "WITH ",
                "collect(rib) AS all_structures,",
                "count(rib) AS total_count",
                f"unwind all_structures[0..$limit] as structure",
                "with total_count, structure,  min(structure.rcsb_id) as next_cursor",
                "return collect(properties(structure)), total_count, next_cursor",
            ]
        )

        query = "\n".join(query_parts)

        print("=======Executing filtered structures query:==========")
        print("\033[96m" + query + "\033[0m")
        print("=====================================================")

        with self.adapter.driver.session() as session:

            def _(tx: Transaction | ManagedTransaction):
                structures, total_count, next_cursor = tx.run(query, params).values()[0]
                return (
                   structures,
                   next_cursor,
                   total_count,
                )

            return session.execute_read(_)

    def get_taxa(self, src_host: typing.Literal["source", "host"]) -> list[int]:
        def _(tx: Transaction | ManagedTransaction):
            return (
                tx.run(
                    """//
                    match (p:PhylogenyNode)-[k:{}]->(s:RibosomeStructure)
                    return collect(distinct(p.ncbi_tax_id))
            """.format(
                        src_host
                    )
                )
                .single()
                .value()
            )

        with self.adapter.driver.session() as session:
            return session.execute_read(_)


dbqueries = Neo4jReader()


#! ---- Sample structures query

# match (rib:RibosomeStructure)
# with rib order by rib.rcsb_id desc

# where
# toLower(rib.citation_title) + toLower(rib.rcsb_id) + toLower(rib.pdbx_keywords_text) + apoc.text.join(rib.citation_rcsb_authors, "")  contains 'complex'
# and (2004 <= rib.citation_year and rib.citation_year <= 2020 or rib.citation_year is null)
# and 0.5 < rib.resolution and rib.resolution < 4.0
# and ALL(x in [ "16SrRNA" ] where x in apoc.coll.flatten(collect{ match (rib)-[]-(p:Polymer) return p.nomenclature }))
# and ANY(tax in [2] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}))

# with collect(rib)[0..20] as rib, count(rib) as total_count
# unwind rib as ribosomes

# optional match (l:Ligand)-[]-(ribosomes)
# with collect(PROPERTIES(l)) as ligands, ribosomes, total_count
# match (rps:Protein)-[]-(ribosomes)
# with collect(PROPERTIES(rps)) as proteins, ligands, ribosomes, total_count
# optional match (rna:RNA)-[]-(ribosomes)
# with collect(PROPERTIES(rna)) as rnas, proteins, ligands, ribosomes, total_count

# with apoc.map.mergeList([{proteins:proteins},{nonpolymeric_ligands:ligands},{rnas:rnas},{other_polymers:[]}]) as rest, ribosomes, total_count
# return collect(apoc.map.merge(ribosomes, rest)),  collect(distinct total_count)[0]


#! ---- Sample polymers query [by_structure]

# match (rib:RibosomeStructure)
# with rib order by rib.rcsb_id desc

# where
# toLower(rib.citation_title) + toLower(rib.rcsb_id) + toLower(rib.pdbx_keywords_text) + apoc.text.join(rib.citation_rcsb_authors, "")  contains 'complex'
# and (2004 <= rib.citation_year and rib.citation_year <= 2020 or rib.citation_year is null)
# and 0.5 < rib.resolution and rib.resolution < 4.0
# and ALL(x in [ "16SrRNA" ] where x in apoc.coll.flatten(collect{ match (rib)-[]-(p:Polymer) return p.nomenclature }))
# and ANY(tax in [2] where tax in apoc.coll.flatten(collect{ match (rib)-[:source]-(p:PhylogenyNode)-[:descendant_of*]-(s:PhylogenyNode) return [p.ncbi_tax_id, s.ncbi_tax_id]}))

# with collect(rib) as rib
# unwind rib as ribosomes

# match (poly:Polymer)-[]-(ribosomes)
# with poly order by poly.nomenclature desc
# with collect(poly)[0..50] as poly, count(poly) as pcount
# unwind poly as polys
# return properties(polys), pcount


#! ---- Sample polymers query [by_polymer_class]

# match (poly:Polymer)-[]-(pc:PolymerClass {class_id:"uL4"})
# with poly order by poly.nomenclature desc
# with collect(poly)[0..50] as poly, count(poly) as pcount
# unwind poly as polys
# return properties(polys), pcount
