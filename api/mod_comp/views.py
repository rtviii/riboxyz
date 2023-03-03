from django.shortcuts import render

# Create your views here.
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from pymol import cmd
from subprocess import Popen, PIPE, STDOUT, run

#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯
@api_view(['GET'])
def hello(request):
    return Response("Hello from comp january")

@api_view(['GET'])
def edit_chain(request):

    params     = dict(request.GET)
    structid   = str.upper(params['structid'][0])

    filehandle = os.path.join(RIBETL_DATA, structid + ".cif")

    cmd.select("resi 1-3")
    cmd.create("{}_test".format(structid), "sele")
    tosave=os.path.join(RIBETL_DATA, f"tempchain{structid}.cif")
    print("Saving to ", tosave)
    cmd.save(tosave)

    return Response("Edited chain successfully")

