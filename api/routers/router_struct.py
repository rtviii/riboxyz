from io import BytesIO
import json
import os
from typing import Optional
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from django.db.models.query import QuerySet
from ninja import FilterSchema,Field
from ninja import Router
from pydantic import BaseModel
from ribctl.etl.ribosome_assets import RibosomeAssets
# from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_auth_asym_id, ranged_align_by_polyclass
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from wsgiref.util import FileWrapper
from ninja import Query
from time import time
structure_router = Router()
TAG              = "Structure"

class StructureFilters(BaseModel): 

      year       : Optional[tuple[int,int]] = None
      resolution : Optional[tuple[int,int]] = None


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

@structure_router.post('/list_structures', response=list[RibosomeStructure], tags=[TAG])
def list_structures(request,filters:StructureFilters):
    import concurrent.futures
    def load_struct(rcsb_id:str):
        return RibosomeAssets(rcsb_id).profile().model_dump_json()

    def read_parallel(rcsb_ids:list[str]):
        with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
            futures = [executor.submit(load_struct, f) for f in rcsb_ids]
            return [fut.result() for fut in futures]
    try:
        struct_ids      = RibosomeAssets.list_all_structs()[:10]
        x               = time()
        struct_profiles = read_parallel(struct_ids[:10])
        y               = time()
        print("opened ", len(struct_profiles))
        print("got filters:", filters)

        return JsonResponse({})

    except Exception as e:
        return HttpResponseServerError("Failed to return structure profiles list :", e)


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

