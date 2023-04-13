import json
import os
import sys
from pathlib import Path
import gql_querystrings

p        = Path(__file__).parents[1]

lsu_path = os.path.join(p,'assets','subunit_map_LSU.json')
ssu_path = os.path.join(p,'assets','subunit_map_SSU.json')

LSU_map = {k: v for k, v in json.load(open(lsu_path, 'r')).items()}
SSU_map = {k: v for k, v in json.load(open(ssu_path, 'r')).items()}


RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
sys.path.append(os.environ.get("PYMOL_PATH"))    