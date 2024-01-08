from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
# from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_auth_asym_id, ranged_align_by_polyclass
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
# from rbxz_bend.application import db_connection

from wsgiref.util import FileWrapper

structure_router = Router()
TAG              = "STRUCTURE"

@structure_router.get('/structure_profile', response=RibosomeStructure, tags=[TAG])
def structure_profile(request,rcsb_id:str):
    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    try:
        with open(RibosomeAssets(rcsb_id)._json_profile_filepath(), 'r') as f:
            return JsonResponse(json.load(f))
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))

#TODO
@structure_router.get('/structure_ptc', response=list[RibosomeStructure], tags=[TAG])
def structure_ptc(request,rcsb_id:str):
    ...
        
#TODO
@structure_router.get('/structure_ligands', response=list[RibosomeStructure], tags=[TAG])
def structure_ligands(request,rcsb_id:str):
    ...
