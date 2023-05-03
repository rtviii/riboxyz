import asyncio
from Bio import AlignIO
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
import os
import sys
from prody import MSA, calcShannonEntropy
from api.rbxz_bend.settings import RIBETL_DATA
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import  Protein, ProteinClass
from api.ribctl.lib.msalib import msa_profiles_dict, msa_profiles_dict_prd, msaclass_extend, msaclass_extend_process_sub

msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()

def seq_H_fit_class(base_class: ProteinClass, base_class_msa:MSA, new_seq:str )->dict[ProteinClass, float]:
    """Calculate entropy difference for a given protein class MSA without and with a new sequence. Used as a measure of fit."""
    extended_class = msaclass_extend(base_class,base_class_msa, new_seq)
    H_original     = sum(calcShannonEntropy(base_class_msa))
    H_extended     = sum(calcShannonEntropy(extended_class))
    H_delta        = H_extended - H_original
    return {base_class: H_delta}

async def seq_H_fit_class_multi(new_seq:str, msa_profiles:dict[ProteinClass, MSA], workers:int=10 )->tuple[ProteinClass,dict[ProteinClass, float]]:

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = []
        for class_name, base_msa in msa_profiles.items():
            fut = executor.submit(seq_H_fit_class, class_name, base_msa, new_seq)
            futures.append(fut)
    wait(futures, return_when=ALL_COMPLETED)

    H_fit:dict[ProteinClass, float] = {}
    for f in futures:
        if f.exception() is not None:
            raise f.exception()
        else:
            H_fit.update(f.result())

    max_fit = min(H_fit, key=H_fit.get)
    return max_fit, H_fit

def fit_chain(chain:Protein,):
    eloop = asyncio.get_event_loop()
    best_fit, fit_dict  = eloop.run_until_complete(seq_H_fit_class_multi(chain.entity_poly_seq_one_letter_code_can, msa_profiles))
    print("Chain with nomenclature {} best fits class {} with delta H = {}".format(chain.nomenclature, best_fit, fit_dict[best_fit]))

def iter_all_profiles_proteins():
    total = {}
    for rcsb_id in os.listdir(RIBETL_DATA):
        if len( rcsb_id ) != 4 or rcsb_id =='6OXI':
            continue
        profile = RibosomeAssets(rcsb_id).profile()
        total[rcsb_id] = {
            'proteins'  : len(profile.proteins),
            'classified': 0
        }
        for protein in profile.proteins:
            if len(protein.nomenclature) != 0:
                total[rcsb_id]['classified'] += 1
    return total

def iter_all_profiles_rna():
    total = {}
    for rcsb_id in os.listdir(RIBETL_DATA):
        if len( rcsb_id ) != 4 or rcsb_id =='6OXI':
            continue
        profile = RibosomeAssets(rcsb_id).profile()
        if profile.rnas == None:
            continue

        total[rcsb_id] = {
            'rna'       : len(profile.rnas),
            'classified': 0
        }
        for rna in profile.rnas:
            if len(rna.nomenclature) != 0:
                total[rcsb_id]['classified'] += 1
    return total

 
rcsb_id:str            = sys.argv[1].upper()
R      :RibosomeAssets = RibosomeAssets(rcsb_id)


prot = R.get_prot_by_nomclass('bL12')
if prot == None: raise LookupError()
msa = msaclass_extend_process_sub('bL12', msa_profiles['bL12'], prot.entity_poly_seq_one_letter_code_can)

for msa_seq in msa:
    print(msa_seq)





# d        = iter_all_profiles_rna()
# all      = d.items()
# percents = [tup[1]['classified']/tup[1]['rna'] for tup in all]

# print("percentage 33", len(list(filter(lambda x: x < 0.33, percents))))
# print("percentage 50", len(list(filter(lambda x: x < 0.5, percents))))
# print("percentage 75", len(list(filter(lambda x: x < 0.75, percents))))
# print("average cov", sum(percents)/len(percents))

# TODO: v0 Coverage  for all structs
# TODO: v1 Coverage  for all structs
# TODO: v0/v1 diff

# The test case for whether (1) this works and (2) this is useful is to run the classification on every structure and track:
# - % of chains where class is assigned for the first time (increased coverage)
# - % of chains where the class assigned conflicts with the existing class (increased accuracy)
# - % of chains where the class assigned is the same as the existing class (no change)
# - all of the above with/without phylogenetic correction (i.e. picking the [10] closest organisms available in a given class )
