# Create your views here.
import sys
from typing import List
from django.http import HttpResponse
from rbxz_bend.neoget import _neoget
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA, CYPHER_EXEC
from rbxz_bend.settings import Neo4jConnection
from neo4j import GraphDatabase, Driver, Session, Transaction, Result, ResultSummary
from subprocess import STDOUT, Popen, PIPE, run
from neo4j.graph import Graph, Node, Relationship, Path
from urllib import parse
import logging
# -⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯


# class Neo4jEntity:
# def commit(*args, **kwargs):


@api_view(['GET'])
def hello(request):
    global connection
    return HttpResponse('hi')


def neo4j_commit_last_db_init():
    cypher = """merge (new:LastUpdated {date:localdatetime()})"""
    with Neo4jConnection.driver.session() as session:
        session.write_transaction(lambda tx: tx.run(
            cypher))


def neo4j_commit_last_update(new_structs: List[str]):
    cypher = """
	     MATCH (n:LastUpdated)
         WITH n
         ORDER BY n.date DESC
         LIMIT 1
         merge (new:LastUpdated {date:localdatetime(), new_structs: $new_structs})-[:previous_update]->(n)
	"""
    with Neo4jConnection.driver.session() as session:
        node_id = session.write_transaction(lambda tx: tx.run(
            cypher, new_structs=new_structs))


# -------------------- SRUCT SINGLE
@api_view(['GET'])
def struct_commit_new(request):
    params = dict(request.GET)
    rcsb_id = str.upper(params['rcsb_id'][0])
    return HttpResponse(neo4j_commit_structure(rcsb_id))

# -------------------- SRUCTS PLURAL


def neo4j_commit_structure(rcsb_id: str):
    rcsb_id = str.upper(rcsb_id)
    proc = Popen(["ribxzcli",  "struct", "obtain",
                 f"{rcsb_id}", "--commit"], env=os.environ.copy(), stdout=PIPE, stderr=PIPE)
    proc.wait()

    out, err = proc.communicate()

    logging.basicConfig(
        filename=f'{RIBETL_DATA}/logs/structs.update.log', filemode='a')
    logger = logging.getLogger(__name__)
    syslog = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s %(app_name)s : %(message)s')

    syslog.setFormatter(formatter)
    logger.setLevel(logging.INFO)
    logger.addHandler(syslog)

    logger = logging.LoggerAdapter(logger, {'app_name': 'ribosome.xyz'})
    update_log = logging.Logger("structs.update")
    update_log.log(logging.INFO, "updated with rcsb_id: " + rcsb_id)

    print(out, err)

    return [out, err]


def neo4j_get_all_struct_ids():
    with GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD)) as driver:
        def transaction_fn(tx: Transaction):
            r = tx.run("match (r:RibosomeStructure) return r.rcsb_id;")
            return r.values()
        with driver.session() as session:
            values = session.read_transaction(transaction_fn)
            return [v[0] for v in values]


def neo4j_diff_w_pdb():
    import requests
    from urllib import parse
    rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query?json="
    params = '''
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
    query = requests.get(rcsb_search_api + parse.quote_plus(params))
    query_response = query.json()
    structs = list(map(lambda x: x['identifier'],
                   query_response['result_set']))

    diff = {
        "rcsb": structs,
        "riboxyz": neo4j_get_all_struct_ids(),
        "diff": list(set(structs) - set(neo4j_get_all_struct_ids()))
    }
    return diff


# ?----------------->>> CACHEABLE

@api_view(['GET'])
def structs_all_ids(request):
    return Response(neo4j_get_all_struct_ids())


# TODO: Centralize logging

@api_view(['GET'])
def structs_sync_with_pdb(request):
    import subprocess
    import shlex
    structs = neo4j_diff_w_pdb()['diff']
    updateloop = "for rcsb_id in {}; do ribxzcli struct obtain \\$rcsb_id --commit >> {}; done".format(
        " ".join(structs),
        f'{RIBETL_DATA}/logs/structs.update.log')
    update_split = shlex.split(updateloop)
    
    subprocess.run(update_split, start_new_session=True)
    neo4j_commit_last_update(structs)
    return HttpResponse(f"{update_split}")


@api_view(['GET'])
def structs_diff_pdb(request):
    return Response(neo4j_diff_w_pdb())
