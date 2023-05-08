import asyncio
import json
import multiprocessing
from concurrent.futures import Future
from concurrent.futures import ALL_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
import os
from pprint import pprint
from prody import MSA, calcShannonEntropy
from api.rbxz_bend.settings import RIBETL_DATA
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import Protein, ProteinClass
from api.ribctl.lib.msalib import msadict_get_meta_info, msa_profiles_dict_prd, msaclass_extend, msaclass_extend_temp
from api.ribctl.taxonomy import filter_by_parent_tax

msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()


def seq_H_fit_class(base_class: ProteinClass, base_class_msa: MSA, target_fasta: str, target_auth_asym_id: str, parent_rcsb_id: str) -> dict[ProteinClass, float]:
    """Calculate entropy difference for a given protein class MSA without and with a new sequence. Used as a measure of fit."""
    extended_class = msaclass_extend_temp(
        base_class, base_class_msa, target_fasta, target_auth_asym_id, parent_rcsb_id)
    H_original = sum(calcShannonEntropy(base_class_msa))
    H_extended = sum(calcShannonEntropy(extended_class))
    H_delta = H_extended - H_original
    # pprint("{} class has {} sequences".format(base_class,base_class_msa.numSequences()))
    return {base_class: H_delta}

async def seq_H_fit_class_multi(chain: Protein, msa_profiles: dict[ProteinClass, MSA], workers=None) -> tuple[ProteinClass, Protein]:
    with ThreadPoolExecutor(max_workers=multiprocessing.cpu_count() if workers is None else workers) as executor:
        futures = []

        for class_name, base_msa in msa_profiles.items():
            fut = executor.submit(seq_H_fit_class, class_name, base_msa,
                                  chain.entity_poly_seq_one_letter_code_can, chain.auth_asym_id, chain.parent_rcsb_id)
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)

    print("Fit chain {}.{}:".format(chain.parent_rcsb_id, chain.auth_asym_id))
    H_fit: dict[ProteinClass, float] = {}
    for f in futures:
        if f.exception() is not None:
            raise f.exception()
        else:
            H_fit.update(f.result())

    max_fit = min(H_fit, key=H_fit.get)
    # pprint(H_fit)
    print("Picked class {}({} sequences) with fit {}".format(max_fit, H_fit[max_fit]))
    return max_fit, chain


def classify_chain(chain: Protein, vvv=False):
    return seq_H_fit_class_multi(chain, msa_profiles)


async def struct_classify_chains(chains: list[Protein], vvv=False):
    return await asyncio.gather(*[seq_H_fit_class_multi(chain, msa_profiles) for chain in chains[:5]]) # * <--------


def process_struct(rcsb_id: str):
    R       = RibosomeAssets(rcsb_id.upper())
    profile = R.profile()
    rps     = profile.proteins

    # if more than one assemlby is present, only use the first one
    if len(profile.assembly_map) > 1:
        first_assembly = []
        for chain_inst in profile.assembly_map[0].polymer_entity_instances:
            first_assembly.append(
                chain_inst.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id)
        rps = list(filter(lambda _: _.auth_asym_id in first_assembly, rps))

    loop = asyncio.get_event_loop()
    result = loop.run_until_complete(struct_classify_chains(rps, vvv=True))

    total = {}
    for assigned_class, chain in result:
        total.update({chain.auth_asym_id: assigned_class})

    return total


def nomv2_duplicates():
    for rcsb_id in os.listdir("./api/ribctl/assets/nomenclaturev2"):
        with open(os.path.join("./api/ribctl/assets/nomenclaturev2", rcsb_id), 'rb') as f:
            prof = json.load(f)

        classes = list(prof.values())
        dup = {x for x in classes if classes.count(x) > 1}

        print(rcsb_id, dup)


def ____generate_nomenclature():
    global counter
    counter = 0
    assets_nomenclaturev2 = os.path.join(
        "/home/rxz/dev/docker_ribxz/api/ribctl/assets/nomenclaturev2")
    loop = asyncio.get_event_loop()

    with ProcessPoolExecutor(max_workers=10) as pool:
        # for struct in filter_by_parent_tax(2):
        for struct in ['3J7Z']:

            def nomv2_path(structid): return os.path.join(
                assets_nomenclaturev2, f"{structid.upper()}.json")

            # if os.path.exists(nomv2_path(struct)):
            #     print("Skipping {}".format(struct))
            #     continue

            if len(struct) != 4 or struct == '6OXI':
                continue
            x = pool.submit(process_struct, struct)

            def save_fut_result(rcsb_id):
                def _(f: Future):
                    global counter
                    counter += 1
                    if f.exception() is not None:
                        print(f.exception())
                        raise LookupError()
                    else:
                        r = f.result()
                        with open(nomv2_path(rcsb_id), 'w') as outfile:
                            json.dump(r, outfile, indent=4)
                return _

            x.add_done_callback(save_fut_result(struct))


if __name__ == "__main__":
    process_struct('3J7Z')


# -------------------------------------------------------------------------------------
# loop.run_until_complete(struct_classify_chains(R.profile().proteins, loop, vvv=True))
# with ProcessPoolExecutor(max_workers=10) as pool:
#     coroutines = []
#     loop.run_until_complete(asyncio.gather(*coroutines))

# =============================

# eloop = asyncio.get_event_loop()
# x = time.time()
# for i in range(10):
#     eloop.run_until_complete(seq_H_fit_class_multi(bl12.entity_poly_seq_one_letter_code_can, msa_profiles))
# y = time.time()
# print("Elapsed: {}".format(y-x))


# d = classify_chains(R.profile().proteins)
# pprint(d)


# def iter_all_profiles_proteins():
#     total = {}
#     for rcsb_id in os.listdir(RIBETL_DATA):
#         if len( rcsb_id ) != 4 or rcsb_id =='6OXI':
#             continue
#         profile = RibosomeAssets(rcsb_id).profile()
#         total[rcsb_id] = {
#             'proteins'  : len(profile.proteins),
#             'classified': 0
#         }
#         for protein in profile.proteins:
#             if len(protein.nomenclature) != 0:
#                 total[rcsb_id]['classified'] += 1
#     return total

    #     if len( rcsb_id ) != 4 or rcsb_id =='6OXI':
    #         continue

    #     profile = RibosomeAssets(rcsb_id).profile()
    #     if profile.rnas == None:
    #         continue

    #     total[rcsb_id] = {
    #         'rna'       : len(profile.rnas),
    #         'classified': 0
    #     }
    #     for rna in profile.rnas:
    #         if len(rna.nomenclature) != 0:
    #             total[rcsb_id]['classified'] += 1
    # return total


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
