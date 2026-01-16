from pprint import pprint
import typing
import sys
sys.path.append('/home/rtviii/dev/riboxyz')
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from ribctl.lib.types.polymer import (
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
from typing import Optional, List, Literal, Tuple, Union

class StructureFilterParams(BaseModel):

    cursor          : Optional[str]                                          = None
    limit           : int                                                    = Field(default=20, ge=1, le=100)
    year            : Optional[tuple[Optional[int], Optional[int]]]          = None
    search          : Optional[str]                                          = None
    resolution      : Optional[tuple[Optional[float], Optional[float]]]      = None
    polymer_classes : Optional[List[PolynucleotideClass | PolypeptideClass]] = None
    source_taxa     : Optional[List[int]]                                    = None
    host_taxa       : Optional[List[int]]                                    = None
    subunit_presence: Optional[Literal["SSU+LSU", "LSU", "SSU"]]             = None

class PolymersFilterParams(BaseModel):

    cursor           : Optional[ Union[Tuple[Optional[str], Optional[str]], List[Optional[str]], str] ] = None
    limit            : int                                                                              = Field(default=20, ge=1, le=100)
    year             : Optional[Tuple[Optional[int], Optional[int]]]                                    = None
    search           : Optional[str]                                                                    = None
    resolution       : Optional[Tuple[Optional[float], Optional[float]]]                                = None
    polymer_classes  : Optional[List[Union[PolynucleotideClass, PolypeptideClass]]]                     = None
    source_taxa      : Optional[List[int]]                                                              = None
    host_taxa        : Optional[List[int]]                                                              = None
    subunit_presence: Optional[Literal["SSU+LSU", "LSU", "SSU"]]                                        = None

    current_polymer_class: Optional[Union[PolynucleotideClass, PolypeptideClass]] = None
    uniprot_id           : Optional[str]                                          = None
    has_motif            : Optional[str]                                          = None

    def get_cursor(self) -> Optional[Tuple[Optional[str], Optional[str]]]:
        if self.cursor is None:
            return None
        if isinstance(self.cursor, str):
            return (self.cursor, self.cursor)
        if isinstance(self.cursor, (list, tuple)) and len(self.cursor) == 2:
            return (self.cursor[0], self.cursor[1])
        raise ValueError("Invalid cursor format")

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

    def ligands_per_structure(self):
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                    result =tx.run(
                        """MATCH (n:Ligand)-[]-(r:RibosomeStructure)
WITH n.chemicalId as chemicalId, collect(r.rcsb_id) as structs
RETURN {chemicalId: chemicalId, structures: structs}"""
                    )
                    return [record[0] for record in result]
            return session.execute_read(_)


    def ligands_in_structure(self, rcsb_id):
        print("GOT HERE", rcsb_id)
        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                result = tx.run(
                    """MATCH (n:Ligand)-[]-(R:RibosomeStructure {rcsb_id: $rcsb_id})
                    RETURN properties(n)
                    """,{
                        "rcsb_id": rcsb_id.upper()
                    }
                )
                k = [record[0] for record in result]
                print("OGT>>>>>>>>>>>>>>>>>",k)
                return k
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

    def list_structs_filtered(self, filters: StructureFilterParams):
        query_parts = ["MATCH (rib:RibosomeStructure)"]
        where_clauses = []
        params = {"limit": filters.limit}

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
                "WITH collect(rib) AS all_structures,",
                "     count(rib) AS total_count",
            ]
        )

        # Apply cursor filter after collecting all structures
        if filters.cursor:
            query_parts.append("WITH all_structures, total_count,")
            query_parts.append(
                "     [s IN all_structures WHERE s.rcsb_id < $cursor] AS filtered_structures"
            )
            params["cursor"] = filters.cursor
        else:
            query_parts.append(
                "WITH all_structures AS filtered_structures, total_count"
            )

        query_parts.extend(
            [
                "WITH filtered_structures[0..{}] AS structures,".format(params["limit"]),
                "     total_count,",
                "     CASE WHEN size(filtered_structures) > {}".format(params["limit"]),
                "          THEN filtered_structures[{}].rcsb_id".format(params["limit"]),
                "          ELSE null",
                "     END AS next_cursor",
                "RETURN",
                "     [struct IN structures | properties(struct)] AS structures,",
                "     total_count,",
                "     next_cursor",
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

    def list_polymers_filtered(self, filters: PolymersFilterParams):
        # First stage: Filter structures
        structure_query_parts = ["MATCH (rib:RibosomeStructure)"]
        structure_where_clauses = []
        params = {"limit": filters.limit}

        if filters.search:
            structure_where_clauses.append(
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
                structure_where_clauses.append(
                    "(datetime(rib.deposition_date).year >= $year_start AND rib.deposition_date IS NOT NULL)"
                )
                params["year_start"] = start
            if end is not None:
                structure_where_clauses.append(
                    "(datetime(rib.deposition_date).year <= $year_end AND rib.deposition_date IS NOT NULL)"
                )
                params["year_end"] = end

        if filters.resolution:
            start, end = filters.resolution
            if start is not None:
                structure_where_clauses.append("rib.resolution > $resolution_start")
                params["resolution_start"] = float(start)
            if end is not None:
                structure_where_clauses.append("rib.resolution < $resolution_end")
                params["resolution_end"] = float(end)

        if filters.polymer_classes:
            structure_where_clauses.append(
                "ALL(x IN $polymer_classes WHERE x IN apoc.coll.flatten(collect{ MATCH (rib)-[]-(p:Polymer) RETURN p.nomenclature }))"
            )
            params["polymer_classes"] = [pc.value for pc in filters.polymer_classes]

        if filters.source_taxa:
            structure_where_clauses.append(
                "EXISTS{ MATCH (rib)-[:belongs_to_lineage_source]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $source_taxa }"
            )
            params["source_taxa"] = filters.source_taxa

        if filters.host_taxa:
            structure_where_clauses.append(
                "EXISTS{ MATCH (rib)-[:belongs_to_lineage_host]-(p:PhylogenyNode) WHERE p.ncbi_tax_id IN $host_taxa }"
            )
            params["host_taxa"] = filters.host_taxa

        if filters.subunit_presence:
            if filters.subunit_presence == "SSU+LSU":
                structure_where_clauses.append(
                    '"lsu" IN rib.subunit_presence AND "ssu" IN rib.subunit_presence'
                )
            elif filters.subunit_presence == "LSU":
                structure_where_clauses.append(
                    '"lsu" IN rib.subunit_presence AND NOT "ssu" IN rib.subunit_presence'
                )
            elif filters.subunit_presence == "SSU":
                structure_where_clauses.append(
                    '"ssu" IN rib.subunit_presence AND NOT "lsu" IN rib.subunit_presence'
                )

        if structure_where_clauses:
            structure_query_parts.append(
                "WHERE " + " AND ".join(structure_where_clauses)
            )

        structure_query_parts.extend(
            [
                "WITH collect(rib) AS filtered_structures",
                "WITH filtered_structures, size(filtered_structures) AS total_structures_count",
            ]
        )

        # Second stage: Match and filter polymers
        query_parts = [
            "UNWIND filtered_structures AS rib",
            "MATCH (rib)-[]-(poly:Polymer)",
            "WHERE poly.assembly_id = 0",
        ]

        polymer_where_clauses = []

        if filters.current_polymer_class:
            polymer_where_clauses.append("$current_polymer_class in poly.nomenclature")
            params["current_polymer_class"] = filters.current_polymer_class.value

        if filters.uniprot_id:
            polymer_where_clauses.append(" $uniprot_id in poly.uniprot_accession")
            params["uniprot_id"] = filters.uniprot_id

        if filters.has_motif:
            polymer_where_clauses.append(
                "poly.entity_poly_seq_one_letter_code_can CONTAINS $has_motif"
            )
            params["has_motif"] = filters.has_motif

        if polymer_where_clauses:
            query_parts.append("AND " + " AND ".join(polymer_where_clauses))

        query_parts.extend(
            [
                "WITH total_structures_count, rib, collect(poly) as rib_polymers",
                "WITH total_structures_count, collect({rib: rib, polymers: rib_polymers}) AS all_data",
                "WITH total_structures_count, all_data, reduce(acc = 0, data IN all_data | acc + size(data.polymers)) AS total_polymers_count",
                "UNWIND all_data AS data",
                "WITH total_structures_count, total_polymers_count, data.rib AS rib, data.polymers AS rib_polymers",
                "UNWIND rib_polymers AS poly",
                "WITH total_structures_count, total_polymers_count, rib, poly",
                "ORDER BY rib.rcsb_id DESC, poly.auth_asym_id ASC",
            ]
        )

        # Apply cursor filter
        cursor = filters.get_cursor()
        if cursor and cursor[0] is not None:
            query_parts.extend(
                [
                    "WITH total_structures_count, total_polymers_count, rib, poly",
                    "WHERE rib.rcsb_id < $cursor_rcsb_id OR (rib.rcsb_id = $cursor_rcsb_id AND poly.auth_asym_id > $cursor_auth_asym_id)",
                ]
            )
            params["cursor_rcsb_id"] = cursor[0]
            params["cursor_auth_asym_id"] = cursor[1] if cursor[1] is not None else ""

        query_parts.extend(
            [
                f"WITH total_structures_count, total_polymers_count, rib, poly",
                f"WITH total_structures_count, total_polymers_count, collect({{rib: rib, poly: poly}}) AS all_polymers",
                f"WITH total_structures_count, total_polymers_count, all_polymers[0..{filters.limit}] AS limited_polymers, size(all_polymers) > {filters.limit} AS has_more",
                "UNWIND limited_polymers AS polymer_data",
                "WITH total_structures_count, total_polymers_count, polymer_data.rib AS rib, polymer_data.poly AS poly, has_more",
                "ORDER BY rib.rcsb_id DESC, poly.auth_asym_id ASC",
                "WITH total_structures_count, total_polymers_count, collect(poly) AS polymers, has_more,",
                "     collect(rib.rcsb_id) AS rcsb_ids, collect(poly.auth_asym_id) AS auth_asym_ids",
                "WITH polymers, total_structures_count, total_polymers_count, has_more,",
                "     rcsb_ids[size(rcsb_ids)-1] AS last_rcsb_id, auth_asym_ids[size(auth_asym_ids)-1] AS last_auth_asym_id",
                "RETURN",
                "    polymers,",
                "    total_structures_count,",
                "    total_polymers_count,",
                "    CASE WHEN has_more THEN {rcsb_id: last_rcsb_id, auth_asym_id: last_auth_asym_id} ELSE null END AS next_cursor",
            ]
        )

        query = "\n".join(structure_query_parts + query_parts)

        print("\n\033[96m" + query + "\033[0m\n")

        with self.adapter.driver.session() as session:
            def _(tx: Transaction | ManagedTransaction):
                result = tx.run(query, params)
                record = result.single()
                if record:
                    polymers, total_structures_count, total_polymers_count, next_cursor = record
                    return (
                        polymers or [],  # Ensure we always return a list, even if empty
                        (
                            (next_cursor["rcsb_id"], next_cursor["auth_asym_id"])
                            if next_cursor
                            else (None, None)
                        ),
                        total_polymers_count or 0,
                        total_structures_count or 0,
                    )
                else:
                    # If no record is returned, provide default values
                    return ([], (None, None), 0, 0)

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
            read = session.execute_read(_)
            print(read)
            return read

dbqueries = Neo4jReader()

