
from typing import Callable
from pydantic import BaseModel
from neo4j import Transaction, ManagedTransaction
from neo4j.graph import Node
from ribctl.lib.libmsa import PhylogenyRank
from ribctl.lib.schema.types_ribosome import PhylogenyNode










# Create constraint on phylogenyNode that ncbi_tax_id is unique
def node__phylogeny(phylogeny_obj:PhylogenyNode)->Callable[[Transaction | ManagedTransaction], Node ]:
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
            merge (phylonode:PhylogenyNode {
                ncbi_tax_id: $ncbi_tax_id,
                scientific_name: $scientific_name,
                rank: $rank,
              })
          return phylonode
        """, **phylogeny_obj.model_dump()).single(strict=True)['phylonode']
    return _

def link__phylogeny(phylogeny_node_neo4j:Node, )->Callable[[Transaction | ManagedTransaction], Node ]:
    
    def _(tx: Transaction | ManagedTransaction):
        return tx.run("""//
            merge (phylonode:PhylogenyNode {
                ncbi_tax_id: $ncbi_tax_id,
                scientific_name: $scientific_name,
                rank: $rank,
              })
          return phylonode
        """, **phylogeny_obj.model_dump()).single(strict=True)['phylonode']
    return _