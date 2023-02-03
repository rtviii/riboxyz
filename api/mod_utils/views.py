# Create your views here.
import sys
from rest_framework.decorators import api_view
from rest_framework.response import Response
import os
from rbxz_bend.settings import  NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA
import subprocess
from subprocess import Popen, PIPE, STDOUT, run
#-⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯⋅⋱⋰⋆⋅⋅⋄⋅⋅∶⋅⋅⋄▫▪▭┈┅✕⋅⋅⋄⋅⋅✕∶⋅⋅⋄⋱⋰⋯⋯⋯


@api_view(['GET'])
def hello(request):
    return Response("Hi again, again.")

@api_view(['GET'])
def struct_assets(request):
    result = subprocess.run([ "ribxzcli", "struct", "show", "7UNW", "--files", "--db"], capture_output=True, text=True)
    print("\nSTDOUT:",result.stdout)
    print("\nSTDERR:",result.stderr)
    return Response(str(result.stderr + result.stdout))

@api_view(['GET'])
def pull_struct_pdb(request):
    params   = dict(request.GET)
    structid = str.upper(params['rcsb_id'][0])
    print("Attempting to acquire PDB for ", structid)

    os.environ["RIBETL_DATA"]         = RIBETL_DATA
    os.environ["EXTRACT_BSITES_PY"]   = "/home/backend/ingress/scripts/extract_bsites.py"
    os.environ["RENDER_THUMBNAIL_PY"] = "/home/backend/ingress/scripts/render_thumbnail.py"
    os.environ["COMMIT_STRUCTURE_SH"] = "/home/backend/ingress/scripts/commit_structure.sh"
    os.environ["SPLIT_RENAME_PY"]     = "/home/backend/ingress/scripts/split_rename.py"

    proc = Popen([ "/home/backend/ingress/src/update_riboxyz.ts",  
                  "--pythonbin", "/opt/venv/bin/python3",
                  "--ingress","commit",
                  "--structure",f"{structid}"], env=os.environ.copy(), stdout=PIPE)

                     

    print([proc.stdout, proc.stdin, proc.stderr])

    sys.stdout.flush()
    return Response([proc.stdout, proc.stdin, proc.stderr])

@api_view(['GET'])
def last_update(request):
    "match (n:Update) return n order by n.date desc  limit 1"
    """MATCH (n:Update)
        WITH n
        ORDER BY n.date DESC
        LIMIT 1
        merge (new:Update {date:date("2023-12-12"), new_structs:['firs','sra']})-[:previous_update]->(n)"""
    return Response()


