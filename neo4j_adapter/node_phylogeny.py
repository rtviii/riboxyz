
from pydantic import BaseModel
from neo4j import Transaction, ManagedTransaction
from ribctl.lib.libmsa import PhylogenyRank










def node__phylogeny(node:PhylogenyNode)->Callable[[Transaction | ManagedTransaction], Node ]:
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