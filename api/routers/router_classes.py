from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
# from ribctl.lib.mod_superimpose import pymol_super, ranged_align_by_auth_asym_id, ranged_align_by_polyclass
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from wsgiref.util import FileWrapper

classification_router = Router()
TAG              = "Polymer Classes"

@classification_router.get('/polynucleotide',  tags=[TAG])
def polynucleotide_class(request,rna_class:PolynucleotideClass):
    params    = dict(request.GET)
    print("Got rna_class", rna_class)
    polyclass = params['rna_class'][0]

    agg = []
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeAssets(rcsb_id).get_chain_by_polymer_class(polyclass)
        except Exception as e:
            print(e)

        if x is not None:
            print("Found {} in {}".format(polyclass,rcsb_id))
            agg.append(x)

    return JsonResponse(json.dumps(agg))

# @classification_router.get('/polypeptide',  tags=[TAG])
# def polypeptide_class(request, _class:PolypeptideClass):
#     params    = dict(request.GET)
#     polyclass = params['_class'][0]

#     response = 'hi'
#     return response
