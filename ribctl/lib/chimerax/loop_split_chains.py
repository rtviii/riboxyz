import os
from chimerax.core.commands import run
from concurrent.futures import ThreadPoolExecutor

RIBETL_DATA = os.environ.get("RIBETL_DATA")

for rcsb_id in os.listdir(RIBETL_DATA):
    try:
        run(session, "ribetl {}".format(rcsb_id))
        run(session, "chainsplitter {}".format(rcsb_id))
        run(session, "close all")
    except Exception as e:
        print(e)