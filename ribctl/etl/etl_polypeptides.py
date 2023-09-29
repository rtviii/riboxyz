import json
import os
import pathlib
import re
import pyhmmer
from pyhmmer import hmmsearch,utils
from pyhmmer.plan7 import HMM
from pyhmmer.easel import  Alphabet,DigitalSequenceBlock
from ribctl.lib.ribosome_types.types_ribosome import LifecycleFactorClass, ProteinClass
from ribctl import ASSETS
from fuzzywuzzy import process, fuzz


LSU_map = {k: v for k, v in json.load(open(ASSETS["subunit_map_lsu"], 'r')).items()}
SSU_map = {k: v for k, v in json.load(open(ASSETS["subunit_map_ssu"], 'r')).items()}

#! Classification / "Type-coercion"

def protein_classify(protein:dict)->list[ProteinClass]:
    banregex = r"/\b([ueb][ls]\d{1,2})\b/gi"
    #     check authors's annotations. if classes are present --> use that.
    finds = re.search(
        banregex, protein['rcsb_polymer_entity']['pdbx_description'])

    if (finds != None):
        firstcap: str = finds[0]
        classname: str = firstcap[0].lower(
        ) + firstcap[1].upper() + firstcap[2:]
        return [classname]

    elif protein['pfams'] == None or protein['pfams'] == []:
        return []
    else:
        pfamids = [pfam['rcsb_pfam_accession'] for pfam in protein['pfams']]
        nomenclature = []
        for pfam_id in pfamids:
            [nomenclature.append(kv[0]) if pfam_id in kv[1]
             ['pfamDomainAccession'] else ... for kv in LSU_map.items()]
            [nomenclature.append(kv[0]) if pfam_id in kv[1]
             ['pfamDomainAccession'] else ... for kv in SSU_map.items()]
        return list(set(nomenclature))

def rna_classify(poly_pdbx_description:str|None):
    rna_reg = {
        "5SrRNA"  : r"\b(?<!\.)(5s)\b",
        "5.8SrRNA": r"\b(?<!\.)(5\.8s)\b",
        "12SrRNA" : r"\b(?<!\.)(12s)\b",
        "16SrRNA" : r"\b(?<!\.)(16s)\b",
        "21SrRNA" : r"\b(?<!\.)(21s)\b",
        "23SrRNA" : r"\b(?<!\.)(23s)\b",
        "25SrRNA" : r"\b(?<!\.)(25s)\b",
        "28SrRNA" : r"\b(?<!\.)(28s)\b",
        "35SrRNA" : r"\b(?<!\.)(35s)\b",
    }

    rnatypes = rna_reg.items()
    for i in rnatypes:
        matches = re.search(i[1], poly_pdbx_description if poly_pdbx_description !=None else '', flags=re.IGNORECASE | re.MULTILINE)
        if matches != None:
            return [i[0]]
    return []

# def factor_classify(description: str) -> LifecycleFactorClass | None:
#     """@description: usually polymer['rcsb_polymer_entity']['pdbx_description'] in PDB"""
#     (match, score) = process.extractOne(description,list_PolymericFactorClass, scorer=fuzz.partial_ratio)
#     return None if score != 100 else match
