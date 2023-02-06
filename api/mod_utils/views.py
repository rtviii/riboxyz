# Create your views here.
import sys
from rbxz_bend.neoget import _neoget
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rbxz_bend.settings import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA, CYPHER_EXEC
from neo4j import GraphDatabase, Driver, Session, Transaction, Result, ResultSummary
from subprocess import Popen, PIPE, STDOUT, run
import requests
from urllib import parse
# -⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯


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


# -------------------- SRUCT SINGLE
@api_view(['GET'])
def struct_commit_new_PDB(request):
    params = dict(request.GET)
    rcsb_id = str.upper(params['rcsb_id'][0])
    proc = Popen(["ribxzcli",  "struct", "obtain",
                 f"{rcsb_id}", "--commit"], env=os.environ.copy(), stdout=PIPE)
    print([proc.stdout, proc.stdin, proc.stderr])
    sys.stdout.flush()
    return Response([proc.stdout, proc.stdin, proc.stderr])

# -------------------- SRUCTS PLURAL


def neo4j_get_all_struct_ids():
    with GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD)) as driver:
        def transaction_fn(tx: Transaction):
            r = tx.run("match (r:RibosomeStructure) return r.rcsb_id;")
            return r.values()
        with driver.session() as session:
            values = session.read_transaction(transaction_fn)
            return [v[0] for v in values]

# ?----------------->>> CACHEABLE
@api_view(['GET'])
def structs_all_ids(request):
    return Response(neo4j_get_all_struct_ids())


@api_view(['GET'])
def structs_diff_pdb(request):
    import requests
    from urllib import parse
    rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    params          = '''
    {
      "query": {
        "type": "group",
        "logical_operator": "and",
        "nodes": [
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "contains_phrase",
              "negation": false,
              "value": "RIBOSOME",
              "attribute": "struct_keywords.pdbx_keywords"
            },
            "node_id": 0
          },
          {
            "type": "terminal",
            "service": "text",
            "parameters": {
              "operator": "greater",
              "negation": false,
              "value": 25,
              "attribute": "rcsb_entry_info.polymer_entity_count_protein"
            },
            "node_id": 1
          }
        ],
        "label": "text"
      },
      "return_type": "entry",
      "request_options": {
        "return_all_hits": true
      }
    }
    '''

    query          = requests.get(rcsb_search_api + parse.quote_plus(params))
    query_response = query.json()
    structs        = list(map(lambda x: x['identifier'], query_response['result_set']))

    diff = {
        "rcsb"   : structs,
        "riboxyz": neo4j_get_all_struct_ids(),
        "diff"   : list(set(structs) - set(neo4j_get_all_struct_ids()))
    }
    return Response(diff)
