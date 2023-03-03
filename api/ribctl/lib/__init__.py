import json
import os
import sys

# open both subunit maps:
LSU_map = {k: v for k, v in json.load(open('ribctl/assets/subunit_map_LSU.json', 'r')).items()}
SSU_map = {k: v for k, v in json.load(open('ribctl/assets/subunit_map_SSU.json', 'r')).items()}

RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
sys.path.append(os.environ.get("PYMOL_PATH")) 