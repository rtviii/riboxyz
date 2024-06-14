import os

from dotenv import load_dotenv

load_dotenv("../.env")

NEO4J_URI       = os.environ.get("NEO4J_URI")
NEO4J_PASSWORD  = os.environ.get("NEO4J_PASSWORD")
NEO4J_USER      = os.environ.get("NEO4J_USER")
NEO4J_CURRENTDB = os.environ.get("NEO4J_CURRENTDB")


print("Gpt neo4j uri", NEO4J_URI)
