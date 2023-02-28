from pprint import pprint
import requests
from urllib.parse import urlencode
from rcsb_api.gql_querystrings import *

gql_structs             = lambda rcsb_id: structure_string.replace("$RCSB_ID", rcsb_id.upper())
gql_polymer_entities    = lambda rcsb_id: polymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())
gql_nonpolymer_entities = lambda rcsb_id: nonpolymer_entities_string.replace("$RCSB_ID", rcsb_id.upper())


def query_rcsb_api(gql_string:str):
	reqstring       = "https://data.rcsb.org/graphql?query={}".format(gql_string)
	try:
		resp  = requests.get(reqstring)
		return resp.json()['data']
	except Exception as e:
		print("Could not land request to RCSB API. {}".format(e))

RCSB_ID  = "3J9M"
# structs  = query_rcsb_api(gql_structs(RCSB_ID))
# nonpolys = query_rcsb_api(gql_nonpolymer_entities(RCSB_ID))
polys    = query_rcsb_api(gql_polymer_entities(RCSB_ID))

rna_only ="""
{
  entry(entry_id: "4UG0 ") {
    polymer_entities {
      entity_poly {
        rcsb_entity_polymer_type
      }
      rcsb_polymer_entity_annotation {
        annotation_id
        assignment_version
        description
        name
        provenance_source
        type
      }
    }
  }
}

"""
polys= query_rcsb_api(rna_only)

for p in polys['entry']['polymer_entities']:
      if p['entity_poly']['rcsb_entity_polymer_type'] == "RNA":	
	      print(p)
# pprint()

# for poly in polys['entry']['polymer_entities']:
# 	if poly["rcsb_polymer_entity_annotation"]!= None:
# 		# print("->>>>",poly["rcsb_polymer_entity_annotation"])
# 		l = list(filter( lambda p: p["type"]=="GO", poly["rcsb_polymer_entity_annotation"]))
# 		for i in l:
# 			print(i)

# 			# print(i['name'])





