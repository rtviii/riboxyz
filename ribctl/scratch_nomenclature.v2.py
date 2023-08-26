import asyncio
import json
import multiprocessing
from concurrent.futures import Future
from concurrent.futures import ALL_COMPLETED, ProcessPoolExecutor, ThreadPoolExecutor, wait
import os
from prody import MSA, calcShannonEntropy
from api.rbxz_bend.settings import RIBETL_DATA
from ribctl.ribosome_assets import RibosomeAssets
from ribctl.lib.types.types_ribosome import Protein, ProteinClass, RibosomeStructure
from ribctl.lib.msalib import msa_dict_get_meta_info, msa_dict, msaclass_extend_temp

CRED      = '\033[91m'
CEND      = '\033[0m'
BACTERIAL = ['4V4I', '4WSD', '4V87', '4WQ1', '4V5D', '1NWX', '4V6Z', '6QNR', '5T7V', '4V9J', '5GAE', '4V4Q', '7UVZ', '7ZQ6', '7QV3', '1VY7', '5LZF', '6GXO', '5MYJ', '4V55', '7OF4', '6HMA', '7UW1', '7RQD', '7K53', '6GSL', '4W2G', '6UO1', '7MSM', '4WRO', '4LEL', '5L3P', '4V48', '7N1P', '4V6N', '4WQU', '4WFN', '7K51', '5O61', '8B7Y', '7P3K', '6ORL', '4V6Y', '6XDQ', '6Q9A', '5EL4', '6S0K', '6WDD', '5U4I', '7JSW', '3JCN', '4V4Z', '4WU1', '7RQA', '6SPG', '1VVJ', '7MSZ', '4V6Q', '4V9Q', '4V72', '4V75', '6GSK', '5KCS', '7UPH', '6V3D', '4V8C', '5WE6', '4V4Y', '6SPD', '7MT2', '6WNT', '5HCP', '4WQY', '4V5H', '1ML5', '5EL6', '4WSM', '6OF6', '7PJX', '6XQD', '5DM6', '5VP2', '8A57', '4V4N', '4V6F', '4V7M', '3JAH', '8B0X', '5KPS', '4IOC', '7JT3', '8AYE', '4V97', '3J7Z', '6WD3', '5LZB', '4V9P', '8EKB', '6QNQ', '6ENJ', '5ADY', '4V5Q', '5J88', '4V42', '1NKW', '4LT8', '5ND9', '5LZW', '5DM7', '4W29', '3J5L', '5DFE', '4V6O', '4V5S', '5M1J', '4V51', '4L71', '1VY5', '7K00', '3J9Y', '5LZD', '7ZQ5', '7UVY', '6OM6', '7A5G', '5A9Z', '7OJ0', '8A63', '6WDK', '6GXN', '8EIU', '7SS9', '7N2V', '6BU8', '7JT2', '5CZP', '4XEJ', '6FKR', '4V95', '6BOH', '6VYR', '6ND6', '3JBU', '3J9Z', '4V4A', '5ZET', '4U1U', '6WD7', '6NDK', '6GXM', '4V4H', '7ZOD', '7UVX', '5AFI', '4V5F', '6Y69', '4P6F', '6V3A', '7JSS', '4YZV', '7OII', '7NHL', '6QDW', '6HA8', '4V8T', '5WFK', '6V3B', '4V5Y', '3JCJ', '4V8U', '5IBB', '6OFX', '4IO9', '5LZZ', '7PAU', '7P6Z', '7MSC', '5AKA', '4V7L', '4V5E', '4WQF', '7N2C', '5HL7', '4V5O', '4W2I', '7P7U', '6CFK', '4V5R', '7YLA', '7RQC', '4Y4O', '3BBX', '4V63', '6N1D', '6HTQ', '4V9O', '3DLL', '5UQ7', '5H5U', '4WZO', '4V54', '6O8Y', '6BUW', '6NSH', '4V8O', '6BY1', '6BZ6', '6VZ3', '7M4W', '7ASM', '5UYL', '5WFS', '1NWY', '6WDA', '4V68', '5LZC', '6WD5', '6B4V', '4WZD', '6X7F', '7AZS', '4V7A', '4V7S', '6S0Z', '4UY8', '3JAI', '7A5F', '4V6T', '4V8J', '7MSH', '4V6D', '5IMQ', '5GAF', '5ZEP', '4V8A', '7KGB', '4YPB', '5VYC', '6TMF', '1YIT', '7P7Q', '7JIL', '7N2U', '5LI0', '4V83', '4V4X', '4TUE', '4WT1', '4V7W', '5LZX', '6WDM', '7SSN', '6WDL', '7ST6', '5GAG', '4V8B', '6YS3', '7MD7', '2ZJR', '7UG7', '4LSK', '6VU3', '4V65', '4V70', '7JQC', '6WOO', '6O97', '4WCE', '4V7J', '4V8H', '6OPE', '6WDC', '6HRM', '7NHN', '6CFJ', '6OGF', '6SPB', '4U26', '7RQB', '7U2J', '5MDV', '4V6G', '6DZI', '4V7I', '2RDO', '7K55', '4V5A', '7Q4K', '4V7X', '6BZ7', '7UNW', '5UYN', '7SSO', '7BV8', '8G61', '5U9F', '4V9I', '5WIS', '4WRA', '6GBZ', '4TUB', '6WD4', '5MDZ', '6ORE', '6OSK', '4V9N', '4V4G', '7NHK', '4V90', '5LZT', '4LFZ', '4V6L', '4U1V', '5WF0', '5GAD', '6WNW', '4V89', '7OIG', '5EL7', '4TUD', '6U48', '4V9B', '6VYQ', '7QGH', '6OT3', '6OTR', '4W2F', '7AQD', '8CVJ', '6VYW', '4CSU', '5J30', '4V7V', '6WD0', '4Y4P', '6Z6K', '5EL5', '5U9G', '4U25', '5KCR', '6VYZ', '7LH5', '6DZP', '6OSQ', '7M5D', '5KPX', '7PJY', '7RYG', '6S13', '1VY6', '5AA0', '6X9Q', '7NWT', '4V6E', '6VZ5', '1W2B', '6WDI', '7PJZ', '7OT5', '6OGI', '8CVL', '5J3C', '7M4V', '7RYF', '7OOD',
             '6I0Y', '7JQB', '5ZLU', '5MDY', '7SFR', '5XYM', '7OIF', '4V84', '6VWL', '6VZ2', '7N30', '6WD9', '6SZS', '4V6V', '5NDK', '8CVK', '7BL4', '4V57', '5APO', '7PJW', '4WOI', '5J8B', '4L47', '6ENF', '4V8F', '6O90', '6XHX', '7P7R', '7O5B', '3CF5', '7S1J', '7ZTA', '6XZA', '6WDE', '6OST', '6O9J', '6WD8', '4U27', '4V8E', '6VYU', '6WD1', '7JSZ', '5KPV', '4TUC', '6OSI', '5J4D', '7QGN', '5E7K', '7RQ8', '6OUO', '2ZJP', '4V79', '6TC3', '4V5K', '4V66', '4WFB', '4V4T', '3PIO', '5IMR', '5V7Q', '4V7K', '4V8D', '8G6Y', '6YSU', '7SA4', '5DOY', '4V7Y', '4U24', '4V5J', '5MGP', '4WWW', '6CFL', '4V78', '8FOM', '7PJT', '4V6A', '6GZQ', '5E81', '6BOK', '7ASO', '4V5P', '5UYM', '4V77', '4U67', '7D6Z', '4V71', '5JVG', '5JVH', '4V64', '7SSL', '6WDH', '4WT8', '4V6S', '4V8X', '4V52', '7U2H', '6VZJ', '7QG8', '5LZY', '7P7S', '7M4Y', '6XDR', '7LV0', '7BL5', '6ENU', '6XHV', '4V8Q', '7NSO', '4Z3S', '5APN', '6XQE', '1XBP', '5J4B', '7S1G', '5UQ8', '6VYS', '5IB8', '4V4V', '6XZ7', '6N9E', '6TBV', '7AZO', '6XZB', '6O9K', '5JTE', '7AQC', '4V6K', '4V9S', '6S12', '5NCO', '6GC0', '1YJW', '4V7T', '7PJU', '6WRU', '7UVW', '6ORD', '7OPE', '4V5B', '4V7Z', '6W6P', '7S1K', '6OG7', '6YSS', '6O8Z', '5GAH', '4V74', '4WF1', '4V49', '6O3M', '5UYP', '4P70', '5W4K', '7QGU', '4V4W', '6H58', '6I7V', '4V5C', '1SM1', '6VYT', '4W4G', '4WR6', '4U20', '7B5K', '6CAE', '6GWT', '7PJS', '5VPO', '6YSI', '5UYK', '5NWY', '4V69', '7S1H', '5J4C', '4IOA', '7BL2', '6XHY', '7UNU', '3JCE', '5LZA', '7OIZ', '7NSQ', '3JAG', '7SSD', '4V85', '6XHW', '6VYY', '2J28', '6WNV', '7RYH', '6DNC', '4V8G', '6VWM', '6WDG', '6VYX', '4V8I', '7M4Z', '4V5G', '4BTS', '6C5L', '7TOS', '4WQR', '5IQR', '6X7K', '6GZZ', '6YEF', '7NWW', '5JC9', '6H4N', '7P7T', '7BL3', '3JCD', '6O8W', '5V8I', '5WE4', '6O8X', '7UNV', '5WDT', '4LNT', '4V47', '6WD2', '3J9W', '8A5I', '7K50', '4V6P', '6GXP', '4V9C', '7S1I', '4V4J', '4V9R', '4YBB', '6OJ2', '7F0D', '6NUO', '7JT1', '7NHM', '6WDF', '7ASP', '4V4P', '4V4R', '7ST2', '6S0X', '5IT8', '7S0S', '7OF2', '6NTA', '4V9H', '7U2I', '4V9D', '8G6X', '7UNR', '5NGM', '6VZ7', '4V5L', '6Q98', '3JA1', '8G6W', '4V7U', '6Q97', '8EKC', '4W2E', '6YSR', '4V67', '6CZR', '5O60', '4V5N', '5WIT', '5OT7', '5J7L', '7QH4', '4WFA', '7UVV', '5FDV', '7OTC', '3J8G', '5J8A', '4V6C', '5LZV', '6WDB', '4V9A', '5JU8', '6V39', '7QQ3', '7LVK', '7BL6', '4V4S', '6C4I', '4V7B', '6UCQ', '4V53', '5FDU', '4V9K', '7RQ9', '5VPP', '4ZSN', '7K52', '7NSP', '8FON', '5DOX', '5UYQ', '4V50', '6ND5', '6VWN', '3JBV', '6HA1', '6OGG', '7QV2', '8BUU', '4V9L', '7PAT', '5ZEB', '5NDJ', '5TCU', '7QGR', '6PJ6', '7MT3', '7RQE', '6GC8', '6Q95', '7N31', '6OXA', '4V9M', '6OF1', '5J91', '6G5I', '7D80', '2ZJQ', '4V5M', '4V6R', '7QV1', '6N9F', '4V73', '4WPO', '6WD6', '1VY4', '6YHS', '6WDJ', '5LZE', '6GZX', '5J5B', '7K54', '8G5Z', '6FXC', '6NWY', '6YST', '5V93', '4V56', '8C8X', '7MT7', '7PJV', '7ST7', '5KPW', '4W2H', '6SPF', '7A5J', '4V7P', '4V8N', '7M4X', '4V76', '5MDW', '4TUA', '3PIP', '6GSJ', '5ND8', '7SSW', '6BZ8']

def seq_H_fit_class(base_class: ProteinClass, base_class_msa: MSA, target_fasta: str, target_auth_asym_id: str, parent_rcsb_id: str,
                     omitgaps: bool = True, vvv: bool = False) -> dict[ProteinClass, float]:
    """Calculate entropy difference for a given protein class MSA without and with a new sequence. Used as a measure of fit."""

    extended_class = msaclass_extend_temp(
        base_class, base_class_msa, target_fasta, target_auth_asym_id, parent_rcsb_id)

    if vvv:
        for seq in base_class_msa:
            print(seq)
        print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
        for seq in extended_class:
            print(seq)
        col_0 = base_class_msa.getArray()[:, 0]
        H_original = calcShannonEntropy(base_class_msa, omitgaps=omitgaps)
        print("\n\n")
        print(len(H_original))

    H_original = sum(calcShannonEntropy(base_class_msa, omitgaps=omitgaps))
    H_extended = sum(calcShannonEntropy(extended_class, omitgaps=omitgaps))
    H_delta = H_extended - H_original
    return {base_class: H_delta}

async def seq_H_fit_class_multi(chain: Protein, msa_profiles_dict: dict[ProteinClass, MSA], workers=None, omitgaps: bool = False, vvv: bool = False) -> tuple[ProteinClass, Protein]:
    with ThreadPoolExecutor(max_workers=min(multiprocessing.cpu_count(), len(msa_profiles_dict)) if workers is None else workers) as executor:
        futures = []
        for class_name, base_msa in msa_profiles_dict.items():
            fut = executor.submit(seq_H_fit_class, class_name, base_msa,
                                  chain.entity_poly_seq_one_letter_code_can, chain.auth_asym_id, chain.parent_rcsb_id, omitgaps, vvv)
            futures.append(fut)

    wait(futures, return_when=ALL_COMPLETED)

    H_fit: dict[ProteinClass, float] = {}
    for f in futures:
        if f.exception() is not None:
            raise f.exception()
        else:
            H_fit.update(f.result())

    max_fit_classname = min(H_fit, key=H_fit.get)

    if vvv:
        meta = msa_dict_get_meta_info(msa_profiles_dict)
        print(f"Fitting chain {chain.parent_rcsb_id}.{chain.auth_asym_id} ({CRED + str(chain.nomenclature) + CEND}) of length L={len(chain.entity_poly_seq_one_letter_code_can)}: ")
        for k, v in H_fit.items():
            if k == max_fit_classname:
                print(
                    f"{CRED}{k}\t<{msa_profiles_dict[k].numSequences()} seqs · {len( msa_profiles_dict[k][0])} AAs>:\t{v}{CEND}")
            else:
                print(
                    f"{k}\t<{msa_profiles_dict[k].numSequences()} seqs · {len( msa_profiles_dict[k][0])} AAs>:\t{v}")
        print("Picked class {}({} sequences) with fit {}".format(
            CRED + str(max_fit_classname) + CEND, meta[max_fit_classname]['nseqs'], H_fit[max_fit_classname]))

    return max_fit_classname, chain

async def struct_classify_chains(chains: list[Protein], vvv: bool = False, omitgaps=False):

    # preload msa profiles by organism
    msa_profiles = {}
    for chain in chains:
        msa_profiles[chain.src_organism_ids[0]] = msa_dict(phylogenetic_correction_taxid=chain.src_organism_ids[0])

    # classify chains
    return await asyncio.gather(*[seq_H_fit_class_multi(chain, msa_profiles[chain.src_organism_ids[0]],  omitgaps=omitgaps,vvv=vvv ) for chain in chains])

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
    result = loop.run_until_complete(struct_classify_chains(rps, omitgaps=False))
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

        # if rcsb_id in BACTERIAL:
        # else:
        #     print("removing ", os.path.join("./api/ribctl/assets/nomenclaturev2", rcsb_id))
        #     os.remove(os.path.join("./api/ribctl/assets/nomenclaturev2", rcsb_id))

def _generate_nomenclature(struct_id_list: list[str], gen_all: bool = False):
    assets_nomenclaturev2 = os.path.join(
        "/home/rxz/dev/docker_ribxz/api/ribctl/assets/nomenclaturev2")
    with ProcessPoolExecutor(max_workers=4) as pool:
        id_list = struct_id_list if not gen_all else os.listdir(RIBETL_DATA)
        # for struct in os.listdir(RIBETL_DATA):
        for struct in id_list:
            def nomv2_path(structid):
                return os.path.join(assets_nomenclaturev2, f"{structid.upper()}.json")

            if os.path.exists(nomv2_path(struct)):
                print("Skipping {}. Exists.".format(nomv2_path(struct)))
                continue

            if len(struct) != 4 or struct == '6OXI':
                continue

            x = pool.submit(process_struct, struct)

            def save_fut_result(rcsb_id):
                def _(f: Future):
                    if f.exception() is not None:
                        print(f.exception())
                        raise LookupError("")
                    else:
                        r = f.result()
                        with open(nomv2_path(rcsb_id), 'w') as outfile:
                            json.dump(r, outfile, indent=4)
                return _

            x.add_done_callback(save_fut_result(struct))


if __name__ == "__main__":

    _generate_nomenclature(BACTERIAL)
    # taxids = {}
    # for rcsb_id in BACTERIAL:
    #     taxids.update({rcsb_id: true_taxid})
    #     with open("true_taxids.json", 'w') as outfile:
    #         json.dump(taxids, outfile, indent=4)

    # nomv2_duplicates()


    # chain, _ = RibosomeAssets('4WRO').get_chain_by_auth_asym_id('4E')
    # chain, _ = RibosomeAssets('4WQU').get_chain_by_auth_asym_id('AL')

    # chain, _ = RibosomeAssets('4V87').get_chain_by_auth_asym_id('A3')
    # chain, _ = RibosomeAssets('4V87').get_chain_by_auth_asym_id('AW')

    # chain, _ = RibosomeAssets('4V87').get_chain_by_auth_asym_id('AZ')
    # organism_taxid = chain.src_organism_ids[0]

    # *ᢹ ritual.rxz[ assets/nomenclaturev2 ]  cat 4W29.json  | grep bS6                                             [master]
    #     "AF": "bS6",
    #     "Be": "bS6",
    #     "Bf": "eL41",
    #     "Bg": "eL41",
    #     "Bh": "eL41",
    # *ᢹ ritual.rxz[ assets/nomenclaturev2 ]  cat 4W29.json  | grep bS6                                             [master]


    # chain, _ = RibosomeAssets('4W29').get_chain_by_auth_asym_id('AF') #  bS6
    # chain, _ = RibosomeAssets('4W29').get_chain_by_auth_asym_id('Be') #! bL12  verify why this is mapped to bs6?

    # chain, _ = RibosomeAssets('4W29').get_chain_by_auth_asym_id('Bf') #  UNK
    # chain, _ = RibosomeAssets('4W29').get_chain_by_auth_asym_id('Bg') #  UNK
    # chain, _ = RibosomeAssets('4W29').get_chain_by_auth_asym_id('Bh') #  UNK



    # ᢹ ritual.rxz[ assets/nomenclaturev2 ]  cat 4V49.json | grep el13                                             [master]
    #     "AR": "eL13",
    #     "BF": "eL13",
    #     "AN": "eL13",
    #     "BW": "eL13",

    # ᢹ ritual.rxz[ assets/nomenclaturev2 ]  cat 4V49.json | grep bl25                                             [master]
    #     "AQ": "bL25",
    #     "BT": "bL25",
    #     "AT": "bL25",



    # msa_profiles = msa_profiles_dict_prd(phylogenetic_correction_taxid=organism_taxid)
    # eloop        = asyncio.get_event_loop()

    # eloop.run_until_complete(seq_H_fit_class_multi(chain, msa_profiles, workers=1, vvv=True, omitgaps=False))


# TODO: Rewrite entropy calculation to account for the alphabet and PDB's shenanigans (4W29.Bf has all "X"s even in canonical)
# TODO: non-profile-based extension of msa classes
# TODO: exit port triangulation

# 3028117
# 2902295

# Fitting chain 4WRO.4E (['uS5']):
# Picked class bL36(69 sequences) with fit 1.5259946410375704

# Fitting chain 4WQU.AL (['uL11']):
# Picked class uL4(162 sequences) with fit 1.8214083307453848


# def old_coverage():
#     all      = d.items()
#     percents = [tup[1]['classified']/tup[1]['rna'] for tup in all]

#     print("percentage 33", len(list(filter(lambda x: x < 0.33, percents))))
#     print("percentage 50", len(list(filter(lambda x: x < 0.5, percents))))
#     print("percentage 75", len(list(filter(lambda x: x < 0.75, percents))))
#     print("average cov", sum(percents)/len(percents))

# TODO: v0 Coverage  for all structs
# TODO: v1 Coverage  for all structs
# TODO: v0/v1 diff

# The test case for whether (1) this works and (2) this is useful is to run the classification on every structure and track:
# - % of chains where class is assigned for the first time (increased coverage)
# - % of chains where the class assigned conflicts with the existing class (increased accuracy)
# - % of chains where the class assigned is the same as the existing class (no change)
# - all of the above with/without phylogenetic correction (i.e. picking the [10] closest organisms available in a given class )
