import json
from pprint import pprint
from django.http import HttpResponse, JsonResponse
from django.shortcuts import render

# Create your views here.
from rest_framework.decorators import api_view
from rest_framework.response import Response

from wsgiref.util import FileWrapper
import os
from subprocess import Popen, PIPE, STDOUT, run

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
@api_view(['GET'])
def hello(request):
    return Response("Hello from comp january")

# @api_view(['GET'])
# def edit_chain(request):

#     params     = dict(request.GET)
#     structid   = str.upper(params['structid'][0])
#     filehandle = os.path.join(RIBETL_DATA, structid + ".cif")
#     cmd.select("resi 1-3")
#     cmd.create("{}_test".format(structid), "sele")
#     tosave=os.path.join(RIBETL_DATA, f"tempchain{structid}.cif")
#     print("Saving to ", tosave)
#     cmd.save(tosave)

#     return Response("Edited chain successfully")


@api_view(['GET'])
def get_chain(request):
    params   = dict(request.GET)
    chainid  = params['auth_asym_id'][0]
    structid = str.upper(params['rcsb_id'][0])
    filename = "{}_STRAND_{}.cif".format(structid, chainid)

    file_handle = os.path.join(RIBETL_DATA,structid,"CHAINS", filename)

    document = open(file_handle, 'rb')
    response = HttpResponse(FileWrapper(document), content_type='chemical/x-mmcif')
    response['Content-Disposition'] = 'attachment; filename="{}_subchain_{}.cif"'.format(structid, chainid)
    return response

@api_view(['GET'])
def get_profile(request):
    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    filename    = "{}.json".format(rcsb_id)
    file_handle = os.path.join(RIBETL_DATA,rcsb_id, filename)

    with open(file_handle, 'rb') as f:
        profile = json.load(f)
        pprint(profile)
        # return HttpResponse(profile, content_type="application/json")

        return JsonResponse(profile)

