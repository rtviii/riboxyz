import os
from chimerax.core.commands import run
from concurrent.futures import ThreadPoolExecutor

from ribctl.etl.etl_assets_ops import Assets

RIBETL_DATA = os.environ.get("RIBETL_DATA")

for rcsb_id in os.listdir('/home/rtviii/dev/RIBETL_DATA'):
    Assets(rcsb_id).paths.chains_dir
    movie_dir = os.path.join(DEST_DIR, "{}.mp4".format(rcsb_id))
    print("Attempting to write ", movie_dir)
    print("Listdir:", os.listdir(DEST_DIR))

    if os.path.join(DEST_DIR, rcsb_id, ".mp4") in os.listdir(DEST_DIR):
        print("Skipping {}. Exists".format(rcsb_id))
        continue
    else:
        run(session,"ribmovie {}".format(rcsb_id))