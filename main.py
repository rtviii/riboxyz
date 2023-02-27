import json
from pprint import pprint
from typing import Any, Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from ribctl.lib.types.types_ribosome import RNA, Ligand, Protein, RibosomeAssets, RibosomeStructure


"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""


def init_driver(uri, username, password):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver(
        "neo4j://localhost:7687", auth=(username, password))
    return api


driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")




# ※ ----------------[ 0.RNA Nodes] 

def link__rna_to_struct(rna: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:

    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
match (rna:RNA) where ID(rna)=$ELEM_ID
match (struct:RibosomeStructure {rcsb_id:$PARENT})
merge (rna)-[rnaof:rna_of]-(struct)
return rna, rnaof, struct""",
                      {"ELEM": int(rna.element_id),
                       "PARENT": parent_rcsb_id}).values('rna', 'rnaof', 'struct')

    return _


def link__rna_to_nomclass(rna: Node) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
MATCH (r:RNA) WHERE ID(r)=$ELEM_ID and r.nomenclature[0] IS NOT NULL with r as rna
MATCH (rna_class:RNAClass {class_id:rna.nomenclature[0]}) 
with rna, rna_class 
merge (rna)-[b:belongs_to]-(rna_class)
return rna, b, rna_class""",
                      {"ELEM": int(rna.element_id)}).values('rna', 'b', 'rna_class')
    return _


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




# ※ ----------------[ 1.Structure Nodes] 

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


# ※ ----------------[ 2.Protein Nodes] 

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


# ※ ----------------[ 3.Ligand Nodes] 

def node__ligand(_ligand:Ligand)->Callable[[Transaction | ManagedTransaction], Node ]:
    L = _ligand.dict()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
 merge (ligand:Ligand {
 	chemicalId          : $chemicalId         ,
 	chemicalName        : $chemicalName       ,
 	formula_weight      : $formula_weight     ,
 	pdbx_description    : $pdbx_description   ,
 	number_of_instances : $number_of_instances
        })
       return ligand
        """, **L).single(strict=True)['ligand']
    return _


def link__ligand_to_struct(prot: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]:
    parent_rcsb_id = parent_rcsb_id.upper()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
  match (ligand:Ligand) where ID(ligand)=$ELEM_ID
  match (struct:RibosomeStructure {rcsb_id: $PARENT})
  merge (ligand)<-[contains:contains]-(struct)
  return struct, ligand, contains
""",
                      {"ELEM_ID": int(prot.element_id),
                       "PARENT": parent_rcsb_id}).values('struct', 'ligand', 'contains')
    return _



#########################################################################################################################################


def rna_init(rna:RNA):
    with driver.session() as s:
        node = s.execute_write(node__rna(rna))
        s.execute_write(link__rna_to_nomclass(node))
        s.execute_write(link__rna_to_struct(node, RNA.parent_rcsb_id))

def protein_init(prot:Protein):
    with driver.session() as s:
        node = s.execute_write(node__protein(prot))
        s.execute_write(link__prot_to_struct(node, prot.parent_rcsb_id))
        x = s.execute_write(link__prot_to_nomclass(node))
        pprint(x)


def ligand_init(lig:Ligand):
    with driver.session() as s:
        node = s.execute_write(node__ligand(lig))
        s.execute_write(link__ligand_to_struct(node, RCSB_ID))


RCSB_ID = "5afi"

rib     = RibosomeStructure.from_json_profile(RCSB_ID)
rna     = rib.rnas[0]
prot    = rib.proteins[0]
lig     = rib.ligands[4]

ligand_init(lig)

# protein_init(prot)