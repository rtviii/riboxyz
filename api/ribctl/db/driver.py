from ast import List
from typing import Dict, LiteralString, Self
from neo4j import Driver, GraphDatabase, Result, Transaction
from pyparsing import Any
from ribctl.db.constraints import NODE_CONSTRAINTS 
from ribctl.db.proteins import node__protein_class
from ribctl.db.rna import node__rna_class
from ribctl.lib.types.types_polymer import list_LSU_Proteins, list_SSU_Proteins, list_RNAClass

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""

# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
# ※ ----------------[ 1.Structure Nodes]
# ※ ----------------[ 2.RNA Nodes]
# ※ ----------------[ 3.Protein Nodes]
# ※ ----------------[ 4.Ligand Nodes] 
# ※ ----------------[ 5.Ingress]

def init_driver(uri:str="neo4j://localhost:7687", username:str="neo4j", password="neo4j"):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver(uri, auth=(username, password))
    return api

driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")

class Neo4jDB():

    driver:Driver 
    def __init__(self, **kwargs)->None:
        self.driver = init_driver(**kwargs)

    def see_constraints(self)->list[dict[str,Any]]:
        with self.driver.session() as s:
            r = s.run("""//
            CALL db.constraints;
                  """)
            return r.data()

    def return_any(self)->list[dict[str,Any]]:
        with self.driver.session() as s:
            return s.read_transaction(lambda tx:tx.run("""//
            match (n) return n limit 10;
            """).data())

    def init_constraints(self)->None:
        with self.driver.session() as s:
            tx:Transaction
            for c in NODE_CONSTRAINTS:
                s.write_transaction(lambda tx: tx.run(c))
                
    def init_protein_classes(self):
        with self.driver.session() as s:
            for protein_class in [*list_LSU_Proteins,*  list_SSU_Proteins]:
                s.execute_write(node__protein_class(protein_class))

    def init_rna_classes(self,):
        with self.driver.session() as s:
            for rna_class in list_RNAClass:
                s.execute_write(node__rna_class(rna_class))




