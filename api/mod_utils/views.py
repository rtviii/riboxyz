# Create your views here.
import sys
from typing import List
from django.http import HttpResponse
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rbxz_bend.neo4j.db import Neo4jConnection
from rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA
from neo4j import GraphDatabase, Driver, Session, Transaction, Result, ResultSummary
from subprocess import STDOUT, Popen, PIPE, run
from neo4j.graph import Graph, Node, Relationship, Path
from urllib import parse
from django.shortcuts import render


# def neo4j_diff_w_pdb():
#     import requests
#     from urllib import parse
#     rcsb_search_api = "https://search.rcsb.org/rcsbsearch/v2/query?json="
#     params = '''
#     {
#       "query": {
#         "type": "group",
#         "logical_operator": "and",
#         "nodes": [
#           {
#             "type": "terminal",
#             "service": "text",
#             "parameters": {
#               "operator": "contains_phrase",
#               "negation": false,
#               "value": "RIBOSOME",
#               "attribute": "struct_keywords.pdbx_keywords"
#             },
#             "node_id": 0
#           },
#           {
#             "type": "terminal",
#             "service": "text",
#             "parameters": {
#               "operator": "greater",
#               "negation": false,
#               "value": 25,
#               "attribute": "rcsb_entry_info.polymer_entity_count_protein"
#             },
#             "node_id": 1
#           }
#         ],
#         "label": "text"
#       },
#       "return_type": "entry",
#       "request_options": {
#         "return_all_hits": true
#       }
#     }
#     '''
#     query = requests.get(rcsb_search_api + parse.quote_plus(params))
#     query_response = query.json()
#     structs = list(map(lambda x: x['identifier'],
#                    query_response['result_set']))

#     diff = {
#         "rcsb": structs,
#         "riboxyz": neo4j_get_all_struct_ids(),
#         "diff": list(set(structs) - set(neo4j_get_all_struct_ids()))
#     }
#     return diff


# ?----------------->>> CACHEABLE