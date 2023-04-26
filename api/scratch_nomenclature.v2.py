import asyncio
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
from pprint import pprint
from prody import MSA, calcShannonEntropy
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import  ProteinClass
from api.ribctl.msa.msalib import msa_profiles_dict_prd, prot_class_msa, prot_class_msa_extend_prd


rcsb_id   :str           = "5AFI"
poly_class:ProteinClass   = "bS6"
R         :RibosomeAssets = RibosomeAssets(rcsb_id)

chain = R.get_chain_by_polymer_class(poly_class)
if chain is None:
    raise LookupError() 

new_seq:str  = chain.entity_poly_seq_one_letter_code_can.replace("\n","")
msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()

# for class_name, msa in list(msa_profiles.items()):
#     class_profile  = prot_class_msa(class_name)
#     extended_class = prot_class_msa_extend_prd( class_name, new_seq)
#     H_extended     = sum(calcShannonEntropy(extended_class))
#     H_original     = sum(calcShannonEntropy(msa))
#     H_delta        = H_extended - H_original

#     print(class_name, "\t:", H_delta)

def seq_asses_protclass_H_fit(base_class: ProteinClass, base_class_msa:MSA, new_seq:str )->dict[ProteinClass, float]:
    """Calculate entropy difference for a given protein class MSA without and with a new sequence. Used as a measure of fit."""
    extended_class = prot_class_msa_extend_prd(base_class,base_class_msa, new_seq)

    H_original     = sum(calcShannonEntropy(base_class_msa))
    H_extended     = sum(calcShannonEntropy(extended_class))
    H_delta        = H_extended - H_original
    # print(class_name, "\t:", H_delta)
    return {base_class: H_delta}

async def compare_against_all_classes(new_seq:str, msa_profiles:dict[ProteinClass, MSA], workers:int=10 )->dict[ProteinClass, float]:
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = []
        for class_name, base_msa in msa_profiles.items():
            fut = executor.submit(seq_asses_protclass_H_fit, class_name, base_msa, new_seq)
            futures.append(fut)
    s = wait(futures, return_when=ALL_COMPLETED)
    fs = [f.result() for f in futures]
    pprint(fs)
    return {}


eloop = asyncio.get_event_loop()
eloop.run_until_complete(compare_against_all_classes(new_seq, msa_profiles))
    