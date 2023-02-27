import json
from pprint import pprint
from typing import Any, Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from ribctl.lib.types.types_ribosome import RNA, Protein, RibosomeAssets, RibosomeStructure


"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""


def init_driver(uri, username, password):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver(
        "neo4j://localhost:7687", auth=(username, password))
    return api


driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")


# ※ -------------------------------= [ 1. Create a Node ]

# # Unit of work
# def get_structures(tx, rcsb_id): # (1)
#     result:Result= tx.run("MATCH (ribosome:RibosomeStructure {rcsb_id: $rcsb_id})-[]-(protein:Protein) RETURN ribosome,protein", rcsb_id=rcsb_id)
#     pprint(result.values('ribosome','protein'))
#     return result

# # Open a Session
# with driver.session() as session:
#     # Run the unit of work within a Read Transaction
#     session.execute_read(get_structures, rcsb_id="5AFI")
#     session.close()

# ※ ----------------[ 2. Neo4j Data Types ]
# Exploring Records

# When accessing a record, either within a loop, list comprehension or within a single record, you can use the [] bracket syntax.

# The following example extracts the p value fro


with driver.session() as s:
    res = s.run(
        "MATCH (lig:Ligand)-[rel]-(rib:RibosomeStructure) RETURN lig,rib,rel limit 5")
    # for record in res.values():
    #     print("\033[91m ------------------------ \033[0m")
    #     pprint(p)

    for record in res:
        # print("\033[91m ------------------------ \033[0m")
        node_lig, node_rib, rel_rel = record["lig"], record["rib"], record['rel']
        # pprint(node_lig.id)
        # pprint(rel_rel.type)
        # pprint(rel_rel.items())
        # pprint(rel_rel.start_node)
        # pprint(rel_rel.end_node)
        # pprint(node_lig.labels)
        # pprint(node_lig.items())


# ※ ----------------[ 2. Neo4j Data Types ]

# _rib = RibosomeStructure(rcsb_id="4UG0")
# _rib = RibosomeAssets("4UG0").json_profile()
rib = RibosomeStructure.from_json_profile("4UG0")


def link__rna_to_struct(rna: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (rna:RNA) where ID(rna)=$ELEM 
match (struct:RibosomeStructure {rcsb_id:"4UG0"})
merge (rna)-[rnaof:rna_of]-(struct)
return rna, rnaof, struct""",
                      {"ELEM": int(rna.element_id),
                       "PARENT": parent_rcsb_id}).values('rna', 'rnaof', 'struct')

    return _


def link__rna_to_nomclass(rna: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
MATCH (r:RNA) WHERE ID(r)=$ELEM and r.nomenclature[0] IS NOT NULL with r as rna
MATCH (rna_class:RNAClass {class_id:rna.nomenclature[0]}) 
with rna, rna_class 
merge (rna)-[b:belongs_to]-(rna_class)
return rna, b, rna_class""",
                      {"ELEM": int(rna.element_id)}).values('rna', 'b', 'rna_class')

    return _


#   with newrna, value
#   match(s:RibosomeStructure {rcsb_id: value.rcsb_id})
#   create (newrna)-[:rna_of]->(s);


def node__rna(_rna: RNA) -> Callable[[Transaction | ManagedTransaction], Node]:
    RNA_dict = _rna.dict()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
merge (rna:RNA {
      asym_ids:     $asym_ids,
      auth_asym_id: $auth_asym_id,

      parent_rcsb_id: $parent_rcsb_id,
      nomenclature  : $nomenclature,
      ligand_like   : $ligand_like,

      src_organism_ids   : $src_organism_ids,
      src_organism_names : $src_organism_names,
      host_organism_ids  : $host_organism_ids,
      host_organism_names: $host_organism_names,


      entity_poly_strand_id               :  $entity_poly_strand_id,
      entity_poly_seq_one_letter_code     :  $entity_poly_seq_one_letter_code,
      entity_poly_seq_one_letter_code_can :  $entity_poly_seq_one_letter_code_can,
      entity_poly_seq_length              :  $entity_poly_seq_length,
      entity_poly_polymer_type            :  $entity_poly_polymer_type,
      entity_poly_entity_type             :  $entity_poly_entity_type
  }) on create set rna.rcsb_pdbx_description = CASE WHEN $rcsb_pdbx_description = null then \"null\" else $rcsb_pdbx_description END

  return rna
""", **RNA_dict).single(strict=True)['rna']

    return _


def node__structure(_rib: RibosomeStructure) -> Callable[[Transaction | ManagedTransaction], Record | None]:

    R = _rib.dict()

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
    merge ( struct:RibosomeStructure{
              rcsb_id               : $rcsb_id,
              expMethod             : $expMethod,
              resolution            : $resolution,
              citation_rcsb_authors : $citation_rcsb_authors,
              citation_title        : $citation_title,

              pdbx_keywords     : $pdbx_keywords     ,
              pdbx_keywords_text: $pdbx_keywords_text,

              src_organism_ids           : $src_organism_ids,
              src_organism_names         : $src_organism_names,

              host_organism_ids           : $host_organism_ids,
             host_organism_names         : $host_organism_names

              })

              on create set
              struct.rcsb_external_ref_id   = CASE WHEN $rcsb_external_ref_id = null then \"null\" else $rcsb_external_ref_id   END,
              struct.rcsb_external_ref_type = CASE WHEN $rcsb_external_ref_type = null then \"null\" else $rcsb_external_ref_type END,
              struct.rcsb_external_ref_link = CASE WHEN $rcsb_external_ref_link = null then \"null\" else $rcsb_external_ref_link END,
              struct.citation_pdbx_doi      = CASE WHEN $citation_pdbx_doi = null then \"null\" else $citation_pdbx_doi END,
              struct.citation_year          = CASE WHEN $citation_year = null then \"null\" else $citation_year END
              
              return struct
        """, **R).single()
    return _


rib = RibosomeStructure.from_json_profile("4UG0")
rna = rib.rnas[0]


with driver.session() as s:
    record = s.execute_write(node__rna(rna))
    print(record.element_id)


def rna_init(RNA:RNA):
    with driver.session() as s:
        node = s.execute_write(node__rna(RNA))
        s.execute_write(link__rna_to_nomclass(node))
        s.execute_write(link__rna_to_struct(node, RNA.parent_rcsb_id))

def protein_init(RNA:RNA):
    with driver.session() as s:
        node = s.execute_write(node__rna(RNA))
        s.execute_write(link__rna_to_nomclass(node))
        s.execute_write(link__rna_to_struct(node, RNA.parent_rcsb_id))
    # pprint(record.get("node_id"))
