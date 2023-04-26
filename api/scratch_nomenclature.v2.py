import asyncio
from concurrent.futures import ALL_COMPLETED, ThreadPoolExecutor, wait
import sys
from prody import MSA, calcShannonEntropy
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from api.ribctl.lib.types.types_ribosome import  Protein, ProteinClass
from api.ribctl.msa.msalib import msa_profiles_dict_prd, prot_class_msa_extend_prd

msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()
def seq_asses_protclass_H_fit(base_class: ProteinClass, base_class_msa:MSA, new_seq:str )->dict[ProteinClass, float]:
    """Calculate entropy difference for a given protein class MSA without and with a new sequence. Used as a measure of fit."""
    extended_class = prot_class_msa_extend_prd(base_class,base_class_msa, new_seq)
    H_original     = sum(calcShannonEntropy(base_class_msa))
    H_extended     = sum(calcShannonEntropy(extended_class))
    H_delta        = H_extended - H_original
    return {base_class: H_delta}

async def compare_against_all_classes(new_seq:str, msa_profiles:dict[ProteinClass, MSA], workers:int=10 )->tuple[ProteinClass,dict[ProteinClass, float]]:

    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = []

        for class_name, base_msa in msa_profiles.items():
            fut = executor.submit(seq_asses_protclass_H_fit, class_name, base_msa, new_seq)
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

def fit_chain(chain:Protein):
    eloop = asyncio.get_event_loop()
    best_fit, fit_dict  = eloop.run_until_complete(compare_against_all_classes(chain.entity_poly_seq_one_letter_code_can, msa_profiles))
    print("Chain with nomenclature {} best fits class {} with delta H = {}".format(chain.nomenclature, best_fit, fit_dict[best_fit]))

rcsb_id:str = sys.argv[1].upper()
R:RibosomeAssets = RibosomeAssets(rcsb_id)

# for c in R.profile().proteins:
#     fit_chain(c)


for cls in list_ProteinClass:
    _ = msa_profiles.get(cls)
    if _ == None:
        print("No MSA for {}".format(cls))


# chain = R.get_chain_by_polymer_class(poly_class)
# if chain is None:
#     raise LookupError() 

