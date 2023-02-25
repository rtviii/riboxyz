import json
import os
from pprint import pprint
from Bio.PDB.Structure import Structure
import itertools
import itertools
import numpy as np
from Bio.PDB import MMCIF2Dict, MMCIFIO, FastMMCIFParser
import gemmi 
flatten = itertools.chain.from_iterable
n1      = np.array
import argparse
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))

def make_nom_dict(profile):
	nomdict = {}
	for i in [*profile['rnas'], *profile['proteins']]:
		nomdict[i['auth_asym_id']]  = i['nomenclature']
	return nomdict

def struct_path(pdbid: str, pftype: str)->str:
    if pftype == 'cif':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}.cif")
    elif pftype == 'json':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}.json") 
    elif pftype == 'modified':
        return os.path.join(RIBETL_DATA, pdbid.upper(), f"{pdbid.upper()}_modified.json") 

def open_structure(pdbid: str, pftype:str):
    pdbid = pdbid.upper()
    if pftype == 'cif':
        cifpath = struct_path(pdbid,'cif')
        try:
            return FastMMCIFParser(QUIET=True).get_structure(pdbid, cifpath)
        except Exception as e:
            return f"\033[93m Parser Error in structure {pdbid} \033[0m : {e}"

    if pftype == 'json':
        with open(struct_path(pdbid, 'json'), 'rb') as _:
            return json.load(_)


structs = {}
for struct_folder in os.listdir(RIBETL_DATA):
    try:
          profile= open_structure(struct_folder, 'json')
          print(make_nom_dict(profile))

          structs = {**structs,**{struct_folder: make_nom_dict(profile)}}
    except Exception as e:
        print(e)
        pass

pprint(structs)
print(len(structs))
with open('nomdicts.json', 'w') as fp:
    json.dump(structs, fp, indent=4)