from django.shortcuts import render

# Create your views here.
from rest_framework.decorators import api_view
from rest_framework.response import Response
from neo4j import Result, GraphDatabase
from rbxz_bend.settings import CYPHER_EXEC
import os


def neoget(cypher_string: str):
    uri = os.getenv('NEO4J_URI')
    user = os.getenv('NEO4J_USER')
    password = os.getenv('NEO4J_PASSWORD')
    current_db = os.getenv('NEO4J_CURRENTDB')
    driver = GraphDatabase.driver(uri, auth=(user, password))

    def parametrized_query(tx, **kwargs):
        result: Result = tx.run(cypher_string, **kwargs)
        return result.value()

    with driver.session() as session:
        return session.read_transaction(parametrized_query)



@api_view(['GET'])
def hello(request):
    return Response("Hello from db")

@api_view(['GET'])
def get_struct(request):
    params = dict(request.GET)
    pdbid = str.upper(params['pdbid'][0])
    cypher = """
    match (n:RibosomeStructure{{rcsb_id:"{pdbid}"}})
    optional match (rr:RNA)-[]-(n)
    with n, collect(rr) as rrna
    optional match (rp:Protein)-[]-(n)
    with n, rrna,  collect(rp) as rps
    optional match (l:Ligand)-[]-(n)
    with n, rrna, rps, collect(l) as ligs
    return {{structure: n, ligands: ligs,rnas: rrna, proteins: rps}}
    """.format_map({"pdbid": pdbid})
    return Response(neoget(cypher))


@api_view(['GET'])
def get_envs(request):
    envs = f"""
    hey:))
    got the following envvars:
    uri        = { os.getenv('NEO4J_URI') }
    user       = { os.getenv('NEO4J_USER') }
    password   = { os.getenv('NEO4J_PASSWORD') }
    current_db = { os.getenv('NEO4J_CURRENTDB') }
          """
    return Response(envs)

@api_view(['GET'])
def check_health(request):
    print(f"""
    got the following envvars:
    uri        = { os.getenv('NEO4J_URI') }
    user       = { os.getenv('NEO4J_USER') }
    password   = { os.getenv('NEO4J_PASSWORD') }
    current_db = { os.getenv('NEO4J_CURRENTDB') }
          """)
    number_of_nodes      = "match (n) return count(n)"
    number_of_structures = "match (n:RibosomeStructure) return count(n)"
    return Response({"struct_number": neoget(number_of_structures), "node_number": neoget(number_of_nodes)})

@api_view(['GET'])
def seed_db(request):
    # execute cypher_execs
    return  

