import json
from pprint import pprint
import typing
from django.http import  JsonResponse, HttpResponseServerError
from ninja import Router
from api.ribxz_api.db_queries import dbqueries
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import  RibosomeStructure
from ribctl.lib.libtax import ncbi

structure_router = Router()
TAG              = "Structure"

@structure_router.get('/profile', response=RibosomeStructure, tags=[TAG],)
def structure_profile(request,rcsb_id:str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    try:
        with open(RibosomeAssets(rcsb_id)._json_profile_filepath(), 'r') as f:
            return JsonResponse(json.load(f))
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))


@structure_router.get('/list_structures', response=list[RibosomeStructure], tags=[TAG])
def list_structures(request):
    structs_response = dbqueries.list_structs()
    structures       = list(map(lambda r: RibosomeStructure.model_validate(r), structs_response))
    return structures



# [
#   {
#     value: '2',
#     title: 'Bacteria',
#     children: [
#       {
#         value: '66',
#         title: 'parent 1-0',
#         children: [
#           {
#             value: '44',
#             title: 'my leaf',
#           },
#           {
#             value: '22',
#             title: 'your leaf',
#           },
#         ],
#       },
#     ],
#   },
#   {
#     value: '4',
#     title: 'Eukarya',
#   },
#   {
#     value: '6',
#     title: 'Prokaroyta',
#   }
# ]


@structure_router.get('/list_source_taxa', response=list, tags=[TAG])
def list_source_taxa(request, src_host:typing.Literal["source", "host"]):
    s = dbqueries.get_taxa(src_host)
   
    return s
