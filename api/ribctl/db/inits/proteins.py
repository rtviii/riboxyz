
from typing import Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.types.types_ribosome import Protein
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass

def node__protein(_prot:Protein)->Callable[[Transaction | ManagedTransaction], Node ]:
    P = _prot.dict()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        merge (rp:Protein {
      asym_ids    : $asym_ids,
      auth_asym_id: $auth_asym_id,
          
      parent_rcsb_id                      : $parent_rcsb_id,
      pfam_comments                       : $pfam_comments,
      pfam_descriptions                   : $pfam_descriptions,
      pfam_accessions                     : $pfam_accessions,

      src_organism_ids   : $src_organism_ids,
      src_organism_names : $src_organism_names,
      host_organism_ids  : $host_organism_ids,
      host_organism_names: $host_organism_names,
      
      ligand_like          : $ligand_like,
      uniprot_accession    : $uniprot_accession,
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


def link__prot_to_struct(prot: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (prot:Protein) where ID(prot)=$ELEM_ID
match (struct:RibosomeStructure {rcsb_id:$PARENT})
merge (prot)-[protof:protein_of]-(struct)
return prot, protof, struct
""",
                      {"ELEM_ID": int(prot.element_id),
                       "PARENT": parent_rcsb_id}).values('prot', 'protof', 'struct')
    return _

def link__prot_to_nomclass(prot: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
   match (prot:Protein) WHERE ID(prot)=$ELEM_ID and prot.nomenclature[0] IS NOT NULL
   merge (prot_class:ProteinClass {class_id:prot.nomenclature[0]})
   merge (prot)-[member:member_of]->(prot_class)
   return prot, member,prot_class
""",
                      {"ELEM_ID": int(prot.element_id)}).values('prot', 'member', 'prot_class')
    return _

def node__protein_class(protein_class:str):
    def _(tx:Transaction | ManagedTransaction):
        return tx.run("""//
            merge (protein_class:NomenclatureClass {class_id:$CLASS_ID})
            return protein_class
        """, {"CLASS_ID":protein_class}).single(strict=True)['protein_class']
    return _

def add_protein(driver:Driver,prot:Protein):
    with driver.session() as s:
        node = s.execute_write(node__protein(prot))
        s.execute_write(link__prot_to_struct(node, prot.parent_rcsb_id))
        s.execute_write(link__prot_to_nomclass(node))