from io import BytesIO
import json
import os
import typing
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
# from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_auth_asym_id, ranged_align_by_polyclass
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
# from rbxz_bend.application import db_connection

from wsgiref.util import FileWrapper

structure_router = Router()
TAG              = "STRUCTURE"
# ? ---------------------- TODOS
# ? ---------------------- TODO

@structure_router.get('/structure_profile', response=list[RibosomeStructure], tags=[TAG])
def structure_profile(request,rcsb_id:str):
    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    filename    = "{}.json".format(rcsb_id)
    file_handle = os.path.join(RIBETL_DATA,rcsb_id, filename)

    try:
        with open(file_handle, 'rb') as f:
            profile = json.load(f)
            return JsonResponse(profile)
    except:
        return HttpResponseServerError("Failed to find structure profile {}".format(rcsb_id))
