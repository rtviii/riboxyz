
from typing import Callable
from neo4j import Driver, ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from api.ribctl.lib.types.types_ribosome import RNA

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
def link__rna_to_nomclass(rna: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
MATCH (r:RNA) WHERE ID(r)=$ELEM_ID and r.nomenclature[0] IS NOT NULL with r as rna
MATCH (rna_class:RNAClass {class_id:rna.nomenclature[0]}) 
with rna, rna_class 
merge (rna)-[b:belongs_to]-(rna_class)
return rna, b, rna_class""",
                      {"ELEM_ID": int(rna.element_id)}).values('rna', 'b', 'rna_class')
    return _

def link__rna_to_struct(rna: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (rna:RNA) where ID(rna)=$ELEM_ID
match (struct:RibosomeStructure {rcsb_id:$PARENT})
merge (rna)-[rnaof:rna_of]-(struct)
return rna, rnaof, struct""",
                      {"ELEM_ID": int(rna.element_id),
                       "PARENT": parent_rcsb_id}).values('rna', 'rnaof', 'struct')
    return _

def node__rna_class(rna_class:str):
    def _(tx:Transaction | ManagedTransaction):
        return tx.run("""//
            merge (rna_class:RNAClass {class_id:$CLASS_ID})
            return rna_class
        """, {"CLASS_ID":rna_class}).single(strict=True)['rna_class']
    return _

def add_rna(driver:Driver,rna:RNA):
    with driver.session() as s:
        node = s.execute_write(node__rna(rna))
        s.execute_write(link__rna_to_nomclass(node))
        s.execute_write(link__rna_to_struct(node, rna.parent_rcsb_id))