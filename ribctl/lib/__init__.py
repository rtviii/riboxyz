import json
import os
import sys
sys.path.append(os.environ.get("PYMOL_PATH")) 

print(f'Invoking __init__.py for {__name__}')

# open both subunit maps:
LSU_map = {k: v for k, v in json.load(open('ribctl/data/subunit_map_LSU.json', 'r')).items()}
SSU_map = {k: v for k, v in json.load(open('ribctl/data/subunit_map_SSU.json', 'r')).items()}
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))