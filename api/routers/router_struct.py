import json
from pprint import pprint
from django.http import  JsonResponse, HttpResponseServerError
from ninja import Router
from api.ribxz_api.db_queries import dbqueries
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import  RibosomeStructure

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

# #TODO
# """map (just stream from emdb), mmcif"""
# @structure_router.get('/mmcif',  tags=[TAG])
# def structure_mmcif(request, rcsb_id:str):
#     params      = dict(request.GET)
#     rcsb_id     = str.upper(params['rcsb_id'][0])
    

#     document = open(RibosomeAssets(rcsb_id)._cif_filepath(), 'rb')
#     response = HttpResponse(FileWrapper(document), content_type='chemical/x-mmcif')
#     response['Content-Disposition'] = 'attachment; filename="{}.cif"'.format(rcsb_id)
#     return response

# @structure_router.get('/ptc', response=list[RibosomeStructure], tags=[TAG])
# def structure_ptc(request,rcsb_id:str):
#     ...
        
# @structure_router.get('/ligands', response=list[RibosomeStructure], tags=[TAG])
# def structure_ligands(request,rcsb_id:str):
#     ...

