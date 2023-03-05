from typing import Self
from neo4j import Driver, GraphDatabase

"""Functions of the form create_node_xxxx return a closure over their target [xxxx] because the neo4j expects a 'unit-of-work'
function that takes a transaction as its argument."""

def init_driver(uri:str="neo4j://localhost:7687", username:str="neo4j", password="neo4j"):
    # TODO: Create an instance of the driver here
    api = GraphDatabase.driver(uri, auth=(username, password))
    return api


driver = init_driver("neo4j://localhost:7687", "neo4j", "neo4j")

class Neo4jDB():

    driver:Driver 

    def __init__(self, uri=None, username=None, password=None)->None:
        self.driver = init_driver(uri, username,password)
        return self

# ※ ----------------[ 0.Database  inits: constraints & nomenclature classes]
# ※ ----------------[ 1.Structure Nodes                                    ]
# ※ ----------------[ 2.RNA       Nodes                                    ]
# ※ ----------------[ 3.Protein   Nodes                                    ]
# ※ ----------------[ 4.Ligand    Nodes                                    ] 
# ※ ----------------[ 4. Ingress]


