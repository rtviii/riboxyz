from neo4j.exceptions import AuthError
from rbxz_bend.db.ribosomexyz import ribosomexyzDB
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER

db_connection = ribosomexyzDB(uri=NEO4J_URI,
                          password=NEO4J_PASSWORD,
                          user=NEO4J_USER
                          )