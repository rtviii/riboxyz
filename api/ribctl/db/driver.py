from ast import List
from typing import Dict, LiteralString, Self
from neo4j import Driver, GraphDatabase, Transaction
from pyparsing import Any
from api.ribctl.db.constraints import NODE_CONSTRAINTS

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""

# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
# ※ ----------------[ 1.Structure Nodes                                    ]
# ※ ----------------[ 2.RNA       Nodes                                    ]
# ※ ----------------[ 3.Protein   Nodes                                    ]
# ※ ----------------[ 4.Ligand    Nodes                                    ] 
# ※ ----------------[ 4. Ingress]

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

    def init_constraints(self)->None:
        with self.driver.session() as s:
            tx:Transaction
            for c in NODE_CONSTRAINTS:
                s.write_transaction(lambda tx: tx.run(c))
                




