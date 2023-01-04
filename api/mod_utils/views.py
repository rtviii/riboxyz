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
def update_db(request):
    # process   = Popen([f"""/home/backend/ingress/src/update_riboxyz.ts  \
    #  --env RIBETL_DATA={RIBETL_DATA}                                    \
    #  --env NEO4J_URI={NEO4J_URI}                                        \
    #  --env NEO4J_PASSWORD={NEO4J_PASSWORD}                              \
    #  --env NEO4J_USER={NEO4J_USER}                                      \
    #  --env NEO4J_CURRENTDB={NEO4J_CURRENTDB}                            \
    #  --env EXTRACT_BSITES_PY={INGRESS_SCRIPTS}/extract_bsites.py        \
    #  --env RENDER_THUMBNAIL_PY={INGRESS_SCRIPTS}/render_thumbnail.py    \
    #  --env SPLIT_RENAME_PY={ INGRESS_SCRIPTS }/split_rename.py          \
    #  --env COMMIT_STRUCTURE_SH={INGRESS_SCRIPTS}/commit_structure.sh    \
    #  --pythonbin {PYTHONBIN}                                            \
    #  --structure {rcsb_id}                                              \
    #  --ingress structure,profile,commit 
    #                  """
    #                  ], stdout=PIPE, stderr=STDOUT)

    # rcsb_id   = "7UNW"
    # PYTHONBIN = "/opt/venv/bin/python3"
    # process   = Popen([f"{INGRESS_EXEC}",
    #  f"--env RIBETL_DATA={RIBETL_DATA}                                    ",
    #  f"--env NEO4J_URI={NEO4J_URI}                                        ",
    #  f"--env NEO4J_PASSWORD={NEO4J_PASSWORD}                              ",
    #  f"--env NEO4J_USER={NEO4J_USER}                                      ",
    # #  f"--env NEO4J_CURRENTDB={NEO4J_CURRENTDB}                            ",

    #  f"--env EXTRACT_BSITES_PY={INGRESS_SCRIPTS}/extract_bsites.py        ",
    #  f"--env RENDER_THUMBNAIL_PY={INGRESS_SCRIPTS}/render_thumbnail.py    ",
    #  f"--env SPLIT_RENAME_PY={INGRESS_SCRIPTS}/split_rename.py          ",
    #  f"--env COMMIT_STRUCTURE_SH={INGRESS_SCRIPTS}/commit_structure.sh    ",

    #  f"--pythonbin {PYTHONBIN}                                            ",
    #  f"--structure {rcsb_id}                                              ",
    #  f"--ingress structure,profile,commit"
    #                  ], stdout=PIPE, stderr=STDOUT)

                     
                #   "--env RIBETL_DATA=/home/backend/api/ribetldata/"
                #   "--env EXTRACT_BSITES_PY=/home/backend/ingress/scripts/extract_bsites.py",
                #   "--env RENDER_THUMBNAIL_PY=/home/backend/ingress/scripts/render_thumbnail.py",
                #   "--env SPLIT_RENAME_PY=/home/backend/ingress/scripts/split_rename.py",
                #   "--env COMMIT_STRUCTURE_SH=/home/backend/ingress/scripts/commit_structure.sh",

    os.environ["RIBETL_DATA"]         = RIBETL_DATA
    os.environ["EXTRACT_BSITES_PY"]   = "/home/backend/ingress/scripts/extract_bsites.py"
    os.environ["RENDER_THUMBNAIL_PY"] = "/home/backend/ingress/scripts/render_thumbnail.py"
    os.environ["COMMIT_STRUCTURE_SH"] = "/home/backend/ingress/scripts/commit_structure.sh"
    os.environ["SPLIT_RENAME_PY"]     = "/home/backend/ingress/scripts/split_rename.py"

    proc = Popen([ "/home/backend/ingress/src/update_riboxyz.ts",  
                  "--pythonbin", "/opt/venv/bin/python3",
                  "--ingress","commit",
                  "--structure","7UNW"], env=os.environ.copy(), stdout=PIPE)

                     

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

@api_view(['GET'])
def pull_struct(request):
    process = Popen(["/home/backend/ingress/hi.ts"], stdout=PIPE, stderr=STDOUT)
    print([process.stdout, process.stdin, process.stderr])
    return Response([process.stdout, process.stdin, process.stderr])

