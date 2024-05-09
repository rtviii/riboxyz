from io import BytesIO
import json
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import FilterSchema,Field
from ninja import Router
from api.ribxz_api.driver import  dbqueries
from neo4j_adapter.adapter import Neo4jAdapter
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import  RibosomeStructure, RibosomeStructureMetadatum
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from wsgiref.util import FileWrapper
from ninja.pagination import paginate
from time import time

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

from ninja.pagination import paginate, PaginationBase
from ninja import Schema

@structure_router.post('/list_structures', response=list[RibosomeStructureMetadatum], tags=[TAG])
def list_structures(request, page:int):
    print(dbqueries.list_structs())

    # def load_metadata(rcsb_id:str):
    #     with open(RibosomeAssets(rcsb_id)._json_profile_filepath(), 'r') as infile:
    #         return RibosomeStructure.model_validate_json(infile.read()).metadata()

    # def read_parallel(rcsb_ids:list[str]):
    #     with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
    #         futures = [executor.submit(load_metadata, f) for f in rcsb_ids]
    #         return [fut.result() for fut in futures]

    # structure_profiles = read_parallel(RibosomeAssets.list_all_structs()[:20])
    return []

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

