from typing import Callable, Literal
from neo4j import  Driver, ManagedTransaction, Record, Transaction
from neo4j.graph import Node, Relationship
from neo4j import ManagedTransaction, Transaction
from ribctl.lib.schema.types_ribosome import  NonpolymericLigand, RibosomeStructure



def struct_exists(rcsb_id:str) -> Callable[[Transaction | ManagedTransaction], bool]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
                MATCH (u:RibosomeStructure {rcsb_id: $rcsb_id})
                return COUNT(u) > 0;
                """, parameters={"rcsb_id":rcsb_id}).single().value()
    return _

# def link__structure_to_phylogeny(rib:RibosomeStructure)->Callable[[Transaction | ManagedTransaction], Node]:
#     rib.src_organism_ids
#     rib.host_organism_ids

def link__structure_to_phylogeny(rib:RibosomeStructure, taxid,relationship:Literal['host_organism', 'source_organism'])->Callable[[Transaction | ManagedTransaction], list[ Node ]]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
        match (struct:RibosomeStructure {rcsb_id: $rcsb_id})
        match (phylo:PhylogenyNode {ncbi_tax_id: $tax_id}) 
        merge (struct)<-[:$relationsip]-(phylo)
        return  struct, struct
        """, {
            "rcsb_id"    : rib.rcsb_id,
            "tax_id"     : taxid,
            "relationsip": relationship
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
                  citation_rcsb_authors : $citation_rcsb_authors,
                  citation_title        : $citation_title,

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
              struct.citation_year          = CASE WHEN $citation_year = null then \"null\" else $citation_year END
              return struct
        """, **R).single()
    return _

# Transaction
def node__ligand(_ligand:NonpolymericLigand)->Callable[[Transaction | ManagedTransaction], Node ]:
    L = _ligand.model_dump()
    def _(tx: Transaction | ManagedTransaction):
     return tx.run("""//
MERGE (ligand:Ligand {chemicalId: 
                $chemicalId })
   ON CREATE SET	ligand.chemicalName        = $chemicalName      ,
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
