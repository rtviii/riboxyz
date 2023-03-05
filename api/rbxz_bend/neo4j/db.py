import os
from typing import Union
from neo4j import GraphDatabase, Driver, Session, Transaction, Result, ResultSummary, Neo4jDriver, BoltDriver


class Neo4jDatabaseConnection():
    driver: Union[Neo4jDriver, BoltDriver]

    def __init__(self, uri, user, password):
        _ = GraphDatabase.driver(uri, auth=(user, password))
        if _ == None:
            Exception("Could not connect to Neo4j")
        else:
            self.driver = _

    def txrun(self, cypher, *args, **kwargs):
        def transaction_fn(tx: Transaction):
            r = tx.run(cypher)
            return r
        with self.driver.session() as session:
            return session.read_transaction(transaction_fn)

    def close(self):
        if self.driver is not None:
            self.driver.close()
        else:
            pass

NEO4J_URI       = os.getenv("NEO4J_URI")
NEO4J_PASSWORD  = os.getenv("NEO4J_PASSWORD")
NEO4J_USER      = os.getenv("NEO4J_USER")
NEO4J_CURRENTDB = os.getenv("NEO4J_CURRENTDB")
Neo4jConnection = Neo4jDatabaseConnection(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)