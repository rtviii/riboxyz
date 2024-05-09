from typing import Callable
from neo4j import Driver, ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import Polymer, Protein


# Transaction
def node__protein(poly:Polymer)->Callable[[Transaction | ManagedTransaction], Node ]:
    P = poly.model_dump()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
    merge (rp:Polymer {
      asym_ids    : $asym_ids,
      auth_asym_id: $auth_asym_id,
      assembly_id: $assembly_id,

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
  rp.rcsb_pdbx_description = CASE WHEN $rcsb_pdbx_description = null then \"null\" else $rcsb_pdbx_description END
  return rp
        """, **P).single(strict=True)['rp']
    return _