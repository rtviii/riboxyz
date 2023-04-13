from api.rbxz_bend.db.ribosomexyz import ribosomexyzDB
from api.rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
import fire



db = ribosomexyzDB(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
fire.Fire(db)