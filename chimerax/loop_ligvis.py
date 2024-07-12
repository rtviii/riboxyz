from chimerax.core.commands import run
# from concurrent.futures import ThreadPoolExecutor
import os
import glob

antibiotic_bsites = "/home/rtviii/dev/riboxyz/antibiotic_bsites/"
json_files = glob.glob(os.path.join(antibiotic_bsites, '*.json'))
for file in json_files:
    chemid_rcsb_id = os.path.basename(file).split(".")[0]
    try:
        run(session, "ligvis {}".format(chemid_rcsb_id))
    except Exception as e:
        print(e)
        

        