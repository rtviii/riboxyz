from typing import Callable, Literal
from neo4j import  Driver, ManagedTransaction, Record, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import  NonpolymericLigand, RibosomeStructure

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

def struct_exists(rcsb_id:str) -> Callable[[Transaction | ManagedTransaction], bool]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
                MATCH (u:RibosomeStructure {rcsb_id: $rcsb_id})
                return COUNT(u) > 0;
                """, parameters={"rcsb_id":rcsb_id}).single().value()
    return _

def link__structure_to_phylogeny(rcsb_id:str, taxid,relationship:Literal['host_organism', 'source_organism'])->Callable[[Transaction | ManagedTransaction], list[ Node ]]:

    def _(tx: Transaction | ManagedTransaction):
        if relationship == 'host_organism':
            return tx.run("""//
            match (struct:RibosomeStructure {rcsb_id: $rcsb_id})
            match (phylo:PhylogenyNode {ncbi_tax_id: $tax_id}) 
            merge (struct)<-[organism:host]-(phylo)
            return  struct, phylo
            """, {
                "rcsb_id"    : rcsb_id,
                "tax_id"     : taxid,
            }).values('struct', 'phylo')
        elif relationship == 'source_organism':
            return tx.run("""//
            match (struct:RibosomeStructure {rcsb_id: $rcsb_id})
            match (phylo:PhylogenyNode {ncbi_tax_id: $tax_id}) 
            merge (struct)<-[organism:source]-(phylo)
            return  struct, phylo
            """, {
                "rcsb_id"    : rcsb_id,
                "tax_id"     : taxid,
            }).values('struct', 'phylo')
    return _

# Transaction
def node__structure(_rib: RibosomeStructure) -> Callable[[Transaction | ManagedTransaction], Record | None]:
    R = _rib.model_dump()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        merge ( struct:RibosomeStructure{
                  rcsb_id               : $rcsb_id,
                  expMethod             : $expMethod,
                  resolution            : $resolution,

                  pdbx_keywords     : $pdbx_keywords,
                  pdbx_keywords_text: $pdbx_keywords_text,

                  src_organism_ids           : $src_organism_ids,
                  src_organism_names         : $src_organism_names,

                  host_organism_ids           : $host_organism_ids,
                  host_organism_names          : $host_organism_names,
                  
                  mitochondrial: $mitochondrial
              })

              on create set
              struct.rcsb_external_ref_id   = CASE WHEN $rcsb_external_ref_id = null then \"null\" else $rcsb_external_ref_id   END,
              struct.rcsb_external_ref_type = CASE WHEN $rcsb_external_ref_type = null then \"null\" else $rcsb_external_ref_type END,
              struct.rcsb_external_ref_link = CASE WHEN $rcsb_external_ref_link = null then \"null\" else $rcsb_external_ref_link END,
              struct.citation_pdbx_doi      = CASE WHEN $citation_pdbx_doi = null then \"null\" else $citation_pdbx_doi END,
              struct.citation_year          = CASE WHEN $citation_year = null then \"null\" else $citation_year END,
              struct.citation_title         = CASE WHEN $citation_title = null then \"null\" else $citation_title END,
              struct.citation_rcsb_authors  = CASE WHEN $citation_rcsb_authors = null then \"null\" else $citation_rcsb_authors END
              return struct
        """, **R).single()
    return _

# Transaction
def node__ligand(_ligand:NonpolymericLigand)->Callable[[Transaction | ManagedTransaction], Node ]:
    L = _ligand.model_dump()
    def _(tx: Transaction | ManagedTransaction):
     return tx.run("""//
    MERGE (ligand:Ligand {chemicalId: $chemicalId })
    ON CREATE SET ligand.chemicalName        = $chemicalName      ,
 	                ligand.formula_weight      = $formula_weight    ,
 	                ligand.pdbx_description    = $pdbx_description  ,
 	                ligand.number_of_instances = $number_of_instances
       RETURN ligand       
        """, **L).single(strict=True)['ligand']
    return _

# Transaction
def link__ligand_to_struct(prot: Node, parent_rcsb_id: str) -> Callable[[Transaction | ManagedTransaction], list[list[Node | Relationship]]]: 
    parent_rcsb_id = parent_rcsb_id.upper()
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
  match (ligand:Ligand) where ELEMENTID(ligand)=$ELEM_ID
  match (struct:RibosomeStructure {rcsb_id: $PARENT})
  merge (ligand)<-[contains:contains]-(struct)
  return struct, ligand, contains
""",
                      {"ELEM_ID": prot.element_id, "PARENT": parent_rcsb_id}).values('struct', 'ligand', 'contains')
    return _

