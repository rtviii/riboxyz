import os
from chimerax.core.commands import run
from concurrent.futures import ThreadPoolExecutor

DEST_DIR = "/home/rtviii/dev/riboxyz/chimerax/movies"

for rcsb_id in ["7K00", "4UG0"]:
    movie_dir = os.path.join(DEST_DIR, "{}.mp4".format(rcsb_id))
    print("Attempting to write ", movie_dir)
    print("Listdir:", os.listdir(DEST_DIR))

    if os.path.join(DEST_DIR, rcsb_id, ".mp4") in os.listdir(DEST_DIR):
        print("Skipping {}. Exists".format(rcsb_id))
        continue
    else:
        run(session,"ribmovie {}".format(rcsb_id))