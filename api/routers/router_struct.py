from io import BytesIO
import json
import os
from typing import Optional
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from django.db.models.query import QuerySet
from pympler.asizeof import asizeof
from ninja import FilterSchema,Field
from ninja import Router
from pydantic import BaseModel
from ribctl.etl.ribosome_assets import RibosomeAssets
# from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_auth_asym_id, ranged_align_by_polyclass
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, RibosomeStructure, RibosomeStructureMetadatum
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
import concurrent.futures
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




class StructureFilters(FilterSchema): 

      year       : Optional[tuple[int,int]] = None
      resolution : Optional[tuple[float,float]] = None
      contains   : Optional[str]            = None


from ninja.pagination import paginate, PaginationBase
from ninja import Schema


class StructurePagination(PaginationBase):
    # only `skip` param, defaults to 5 per page
    class Input(Schema):
        skip: int


    class Output(Schema):
        items: list[RibosomeStructure] # `items` is a default attribute
        total: int
        per_page: int = 5

    def paginate_queryset(self, queryset, pagination: Input, **params):
        skip = pagination.skip
        print("Got queryset", queryset.count)
        return {
            'items'   : queryset[skip : skip + 5],
            'total'   : queryset.count(),
            'per_page': 5,
        }




@structure_router.post('/list_structures', response=list[RibosomeStructureMetadatum], tags=[TAG])
@paginate(pagination=StructurePagination)
def list_structures(request):
    def load_struct(rcsb_id:str):
        return RibosomeAssets(rcsb_id).profile().metadata()

    def read_parallel(rcsb_ids:list[str]):
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(load_struct, f) for f in rcsb_ids]
            return [fut.result() for fut in futures]

    # try:
    # struct_ids      = RibosomeAssets.list_all_structs()[:10]
    # struct_profiles = [RibosomeAssets(rcsb_id).profile() for rcsb_id in struct_ids]

    structure_profiles = read_parallel(RibosomeAssets.list_all_structs()[:10])
    return structure_profiles

    # except Exception as e:
    #     return HttpResponseServerError("Failed to return structure profiles list :", e)


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

