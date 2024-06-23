from pprint import pprint
from typing import Callable, Literal
from neo4j import Driver, ManagedTransaction, Record, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import NonpolymericLigand, RibosomeStructure

# Get superkingdom given an rcsb_id
"""
match (n:RibosomeStructure)-[]-(p:PhylogenyNode {ncbi_tax_id:83333}) 
with n, p
match (p)-[:descendant_of*]-(s {rank:"superkingdom"})
return n,s;

"""
# get all structs for a given superkingdom
"""
match (p:PhylogenyNode)-[:descendant_of*]-(s {ncbi_tax_id:2,rank:"superkingdom"})
with p,s
match (p)-[]-(str:RibosomeStructure)
return p.ncbi_tax_id, str.rcsb_id
"""


def struct_exists(rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], bool]:
    def _(tx: Transaction | ManagedTransaction):
        return (
            tx.run(
                """//
                MATCH (u:RibosomeStructure {rcsb_id: $rcsb_id})
                return COUNT(u) > 0;
                """,
                parameters={"rcsb_id": rcsb_id},
            )
            .single()
            .value()
        )

    return _


def link__structure_to_lineage_member(
    rcsb_id: str,
    taxid,
    relationship: Literal["belongs_to_lineage_host", "belongs_to_lineage_source"],
) -> Callable[[Transaction | ManagedTransaction], list[Node]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
            match (struct:RibosomeStructure {{rcsb_id: $rcsb_id}})
            match (phylo:PhylogenyNode {{ncbi_tax_id: $tax_id}}) 
            merge (struct)<-[organism: {} ]-(phylo)
            return  struct, phylo
            """.format(relationship),
            {
                "rcsb_id": rcsb_id,
                "tax_id": taxid,
            },
        ).values("struct", "phylo")

    return _


def link__structure_to_organism(
    rcsb_id: str, taxid, relationship: Literal["host", "source"]
) -> Callable[[Transaction | ManagedTransaction], list[Node]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
            match (struct:RibosomeStructure {{rcsb_id: $rcsb_id}})
            match (phylo:PhylogenyNode {{ncbi_tax_id: $tax_id}}) 
            merge (struct)<-[organism: {} ]-(phylo)
            return  struct, phylo
            """.format(relationship),
            {
                "rcsb_id": rcsb_id,
                "tax_id": taxid,
            },
        ).values("struct", "phylo")

    return _


# Transaction
def node__structure(
    _rib: RibosomeStructure,
) -> Callable[[Transaction | ManagedTransaction], Record | None]:
    R = _rib.model_dump()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
        MERGE (struct:RibosomeStructure{ rcsb_id: $rcsb_id})
        SET
            struct.expMethod = $expMethod,
            struct.resolution = $resolution,
            struct.pdbx_keywords = $pdbx_keywords,
            struct.pdbx_keywords_text = $pdbx_keywords_text,
            struct.src_organism_ids = $src_organism_ids,
            struct.src_organism_names = $src_organism_names,
            struct.host_organism_ids = $host_organism_ids,
            struct.host_organism_names = $host_organism_names,
            struct.mitochondrial = $mitochondrial,
            struct.rcsb_external_ref_id = CASE WHEN $rcsb_external_ref_id IS NULL THEN "null" ELSE $rcsb_external_ref_id END,
            struct.rcsb_external_ref_type = CASE WHEN $rcsb_external_ref_type IS NULL THEN "null" ELSE $rcsb_external_ref_type END,
            struct.rcsb_external_ref_link = CASE WHEN $rcsb_external_ref_link IS NULL THEN "null" ELSE $rcsb_external_ref_link END,
            struct.citation_pdbx_doi = CASE WHEN $citation_pdbx_doi IS NULL THEN "null" ELSE $citation_pdbx_doi END,
            struct.citation_year = CASE WHEN $citation_year IS NULL THEN "null" ELSE $citation_year END,
            struct.citation_title = CASE WHEN $citation_title IS NULL THEN "null" ELSE $citation_title END,
            struct.citation_rcsb_authors = CASE WHEN $citation_rcsb_authors IS NULL THEN "null" ELSE $citation_rcsb_authors END,
            struct.subunit_presence = $subunit_presence
        RETURN struct
               """ , **R).single()

    return _


# Transaction
def node__ligand(
    _ligand: NonpolymericLigand,
) -> Callable[[Transaction | ManagedTransaction], Node]:
    L = _ligand.model_dump()

    def _(tx: Transaction | ManagedTransaction):
        drugbank_id = None
        drugbank_description = None
        try:
            drugbank_id = (
                _ligand.nonpolymer_comp.drugbank.drugbank_container_identifiers.drugbank_id
            )
        except:
            ...
        try:
            drugbank_description = (
                _ligand.nonpolymer_comp.drugbank.drugbank_info.description
            )
        except:
            ...

        return tx.run(
            """//
    MERGE (ligand:Ligand {chemicalId: $chemicalId })
    ON CREATE SET   ligand.chemicalName        = $chemicalName      ,
 	                ligand.formula_weight      = $formula_weight    ,
 	                ligand.pdbx_description    = $pdbx_description  ,
 	                ligand.number_of_instances = $number_of_instances,
                    ligand.drugbank_id          =  CASE WHEN $drugbank_id = null then \"null\" else $drugbank_id END,
                    ligand.drugbank_description =  CASE WHEN $drugbank_description = null then \"null\" else $drugbank_description END
               RETURN ligand       
        """,
            {
                "chemicalId": _ligand.chemicalId,
                "chemicalName": _ligand.chemicalName,
                "formula_weight": _ligand.formula_weight,
                "pdbx_description": _ligand.pdbx_description,
                "number_of_instances": _ligand.number_of_instances,
                "drugbank_id": drugbank_id,
                "drugbank_description": drugbank_description,
            },
        ).single(strict=True)["ligand"]

    return _


# Transaction
def link__ligand_to_struct(
    prot: Node, parent_rcsb_id: str
) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    parent_rcsb_id = parent_rcsb_id.upper()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run(
            """//
  match (ligand:Ligand) where ELEMENTID(ligand)=$ELEM_ID
  match (struct:RibosomeStructure {rcsb_id: $PARENT})
  merge (ligand)<-[contains:contains]-(struct)
  return struct, ligand, contains
""",
            {"ELEM_ID": prot.element_id, "PARENT": parent_rcsb_id},
        ).values("struct", "ligand", "contains")

    return _
