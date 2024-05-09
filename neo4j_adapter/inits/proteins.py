from typing import Callable
from neo4j import Driver, ManagedTransaction, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import Protein

def node__protein(_prot:Protein)->Callable[[Transaction | ManagedTransaction], Node ]:
    P = _prot.model_dump()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        merge (rp:Protein {
      asym_ids    : $asym_ids,
      auth_asym_id: $auth_asym_id,
      assembly_id: $assembly_id,

          
      parent_rcsb_id                      : $parent_rcsb_id,
      pfam_comments                       : $pfam_comments,
      pfam_descriptions                   : $pfam_descriptions,
      pfam_accessions                     : $pfam_accessions,

      src_organism_ids   : $src_organism_ids,
      src_organism_names : $src_organism_names,
      host_organism_ids  : $host_organism_ids,
      host_organism_names: $host_organism_names,
      
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

    print("Got protein node :", prot._element_id)
    print("Got protein node :", prot.element_id)
    print("Attemptint ot link to parent structure with rcsb_id:", parent_rcsb_id)
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (prot:Protein) where ELEMENTID(prot)=$ELEM_ID
match (struct:RibosomeStructure {rcsb_id:$PARENT})
merge (prot)-[protof:protein_of]-(struct)
return prot, protof, struct
""",
                      {"ELEM_ID": prot.element_id,
                       "PARENT": parent_rcsb_id}).values('prot', 'protof', 'struct')
    return _

def link__prot_to_polymer_class(prot: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    print("Attempting to link to polymer class")
    print("Got node id", prot.id)
    print("Got node element_id", prot.element_id)
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
   match (prot:Protein) WHERE ELEMENTID(prot)=$ELEM_ID and prot.nomenclature[0] IS NOT NULL
   match (polymer_class:PolymerClass {class_id:prot.nomenclature[0]})
   merge (prot)-[member:member_of]->(polymer_class)
   return prot, member,polymer_class
""",
                      {"ELEM_ID": prot.element_id}).values('prot', 'member', 'polymer_class')
    return _

def node__polymer_class(polymer_class:str):
    def _(tx:Transaction | ManagedTransaction):
        return tx.run("""//
            merge (polymer_class:PolymerClass {class_id:$CLASS_ID})
            return polymer_class
        """, {"CLASS_ID":polymer_class}).single(strict=True)['polymer_class']
    return _

def add_protein(driver:Driver,prot:Protein):
    with driver.session() as s:
        node = s.execute_write(node__protein(prot))
        s.execute_write(link__prot_to_struct(node, prot.parent_rcsb_id))
        s.execute_write(link__prot_to_polymer_class(node))