from neo4j.exceptions import AuthError
from ribctl.db.ribosomexyz import riboxyzDB
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER

# global db_connection

# try:
#     print("trying -----------")
#     db_connection    = riboxyzDB(uri=NEO4J_URI,password=NEO4J_PASSWORD,user=NEO4J_USER)

# except AuthError as e:

#     print("Could not connect to Neo4j database. Trying to change the default password.")

#     _temp_db_connection= riboxyzDB(uri=NEO4J_URI, password="neo4j",user="neo4j")
#     with _temp_db_connection.driver.session(database='system') as system_s:
#         r = system_s.run("""ALTER CURRENT USER SET PASSWORD FROM "neo4j" TO "ribosomexyz";""")
#         print("Changed the default Neo4j password.")

#     print("trying again -----------")
#     db_connection= riboxyzDB(uri= NEO4J_URI,password = NEO4J_PASSWORD,user= NEO4J_USER)
db_connection= riboxyzDB(uri= NEO4J_URI,password = NEO4J_PASSWORD,user= NEO4J_USER)
# db_connection= riboxyzDB(uri=NEO4J_URI, password="neo4j",user="neo4j")