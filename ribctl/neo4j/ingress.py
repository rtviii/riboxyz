from typing import Any, Callable
from neo4j import GraphDatabase, Driver, ManagedTransaction, Record, Result, Transaction
from neo4j.graph import Node, Relationship
from ribctl.lib.types.types_ribosome import RNA, Ligand, Protein, ProteinClass, RibosomeAssets, RibosomeStructure
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass




def init_driver(NEO4J_URI:str="neo4j://localhost:7687", NEO4J_USERNAME:str="neo4j", NEO4J_PASSWORD:str="neo4j"):
    return GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USERNAME, NEO4J_PASSWORD))


# TODO: write this into a cli, encrypted WHEN/IF people actually need access. Keep auth disabled for now
# def change_default_password(new_password:str): #     with driver.session(database="system") as s:
#         print(s.run("show default database").single())
#     s.run("ALTER CURRENT USER SET PASSWORD FROM $CURRENT_PASSWORD TO $NEW_PASSWORD", {
#         "CURRENT_PASSWORD": "neo4j",
#         "NEW_PASSWORD": "55288",
#     } ).single()
    

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""
### DATABASE INITIALIZATION
def node__protein_class(protein_class:str):
    def _(tx:Transaction | ManagedTransaction):
        return tx.run("""//
            merge (protein_class:NomenclatureClass {class_id:$CLASS_ID})
            return protein_class
        """, {"CLASS_ID":protein_class}).single(strict=True)['protein_class']
    return _

def node__rna_class(rna_class:str):
    def _(tx:Transaction | ManagedTransaction):
        return tx.run("""//
            merge (rna_class:RNAClass {class_id:$CLASS_ID})
            return rna_class
        """, {"CLASS_ID":rna_class}).single(strict=True)['rna_class']
    return _

def init_protein_classes(driver:Driver=init_driver()):
    with driver.session() as s:
        for protein_class in [*list_LSU_Proteins,*  list_SSU_Proteins]:
            s.execute_write(node__protein_class(protein_class))

def init_rna_classes(driver:Driver=init_driver()):
    with driver.session() as s:
        for rna_class in list_RNAClass:
            s.execute_write(node__rna_class(rna_class))

def create_constraints(tx: Transaction | ManagedTransaction):
    tx.run("""//
    CREATE CONSTRAINT IF NOT EXISTS ON (ipro:InterProFamily) ASSERT ipro.family_id  IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (go:GOClass) ASSERT go.class_id IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (q:RibosomeStructure) Assert q.rcsb_id IS UNIQUE;
    CREATE CONSTRAINT IF NOT EXISTS ON (pf:PFAMFamily) assert pf.family_id  is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (lig:Ligand) assert lig.chemicalId is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;
    CREATE CONSTRAINT IF NOT EXISTS ON (nc:NomenclatureClass) assert nc.class_id is unique;
    """)

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

# ※ ----------------[ 4. Ingress]
def commit_init(rna:RNA, driver:Driver=init_driver()):
    with driver.session() as s:
        node = s.execute_write(node__rna(rna))
        s.execute_write(link__rna_to_nomclass(node))
        s.execute_write(link__rna_to_struct(node, RNA.parent_rcsb_id))

def commit_protein(prot:Protein, driver:Driver=init_driver()):
    with driver.session() as s:
        node = s.execute_write(node__protein(prot))
        s.execute_write(link__prot_to_struct(node, prot.parent_rcsb_id))
        x = s.execute_write(link__prot_to_nomclass(node))
        pprint(x)

def commit_ligand(lig:Ligand, driver:Driver=init_driver()):
    with driver.session() as s:
        node = s.execute_write(node__ligand(lig))
        s.execute_write(link__ligand_to_struct(node, RCSB_ID))



