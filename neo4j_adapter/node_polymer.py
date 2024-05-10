from typing import Callable
from neo4j import Driver, ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import RNA, Polymer, Protein 


def node__polymer(poly:Polymer)->Callable[[Transaction | ManagedTransaction], Node ]:
    P = poly.model_dump()
    print("Creating poly node", P)
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
    merge (poly:Polymer {
      assembly_id: $assembly_id,

      asym_ids    : $asym_ids,
      auth_asym_id: $auth_asym_id,

      parent_rcsb_id                      : $parent_rcsb_id,

      src_organism_ids   : $src_organism_ids,
      src_organism_names : $src_organism_names,

      host_organism_ids  : $host_organism_ids,
      host_organism_names: $host_organism_names,
      
      rcsb_pdbx_description: $rcsb_pdbx_description,

      entity_poly_strand_id              : $entity_poly_strand_id,
      entity_poly_seq_one_letter_code    : $entity_poly_seq_one_letter_code,
      entity_poly_seq_one_letter_code_can: $entity_poly_seq_one_letter_code_can,
      entity_poly_seq_length             : $entity_poly_seq_length,
      entity_poly_polymer_type           : $entity_poly_polymer_type,
      entity_poly_entity_type            : $entity_poly_entity_type,

      nomenclature                       : $nomenclature
  })
  on create set
  poly.rcsb_pdbx_description = CASE WHEN $rcsb_pdbx_description = null then \"null\" else $rcsb_pdbx_description END
  return poly
        """, **P).single(strict=True)['poly']
    return _

def link__polymer_to_structure(polymer: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (polymer:Polymer) where ELEMENTID(polymer)=$ELEM_ID
match (struct:RibosomeStructure {rcsb_id:$PARENT})
merge (polymer)<-[poly:Subcomponent]-(struct)
return polymer, poly, struct
""",
                      {"ELEM_ID": polymer.element_id,
                       "PARENT": parent_rcsb_id}).values('polymer', 'poly', 'struct')
    return _

def link__polymer_to_polymer_class(poly: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
   match (poly:Polymer) WHERE ELEMENTID(poly)=$ELEM_ID and poly.nomenclature[0] IS NOT NULL
   match (polymer_class:PolymerClass {class_id:poly.nomenclature[0]})
   merge (poly)-[member:member_of]->(polymer_class)
   return poly, member,polymer_class
    """,
                      {"ELEM_ID": poly.element_id}).values('poly', 'member', 'polymer_class')
    return _

def upsert_polymer_to_protein( polymer_node:Node, protein:Protein):
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        match (n) where elementid(n)=$ELEM_ID 
        with n, $PROT as map
        with n, apoc.map.removeKeys(map, keys(n)) as onlyNewPropsMap
            set n+=onlyNewPropsMap, n:Protein
        return n
    """, {"ELEM_ID": polymer_node.element_id, "PROT": protein.model_dump()}).single()

    return _

def upsert_polymer_to_rna( polymer_node:Node, rna:RNA):
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        match (n) where elementid(n)=$ELEM_ID 
        set n:RNA
        return n
    """, {"ELEM_ID": polymer_node.element_id, "RNA": rna.model_dump()}).single()

    return _