# Create your views here.
import sys
from rbxz_bend.neoget import _neoget
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rbxz_bend.settings import  NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA, CYPHER_EXEC
from neo4j import GraphDatabase,Driver,Session, Transaction, Result, ResultSummary
from subprocess import Popen, PIPE, STDOUT, run
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯


@api_view(['GET'])
def hello(request):
    return Response("Hi again, again.")


@api_view(['GET'])
def last_update(request):
    "match (n:Update) return n order by n.date desc  limit 1"
    """MATCH (n:Update)
        WITH n
        ORDER BY n.date DESC
        LIMIT 1
        merge (new:Update {date:date("2023-12-12"), new_structs:['firs','sra']})-[:previous_update]->(n)"""
    return Response()


#-------------------- SRUCT SINGLE
@api_view(['GET'])
def struct_commit_new_PDB(request):
    params  = dict(request.GET)
    rcsb_id = str.upper(params['rcsb_id'][0])
    proc    = Popen([ "ribxzcli",  "struct", "obtain", f"{rcsb_id}", "--commit"], env=os.environ.copy(), stdout=PIPE)
    print([proc.stdout, proc.stdin, proc.stderr])
    sys.stdout.flush()
    return Response([proc.stdout, proc.stdin, pro.stderr])

#-------------------- SRUCTS PLURAL
@api_view(['GET'])
def structs_get_ids(request):
    with GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD)) as driver:
        def transaction_fn(tx:Transaction):
             r = tx.run("match (r:RibosomeStructure) return r.rcsb_id;")
             return r.values()
        with driver.session() as session:
            values = session.read_transaction(transaction_fn)
            return Response([v[0] for v in values])
