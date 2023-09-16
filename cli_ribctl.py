# TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.
# export DJANGO_SETTINGS_MODULE=api.rbxz_bend                                                                                         [cli]
# import argparse
import argparse
import asyncio
import json
import os
from ribctl import RIBETL_DATA
from ribctl.etl.etl_pipeline import ReannotationPipeline, query_rcsb_api, rcsb_single_structure_graphql
from ribctl.ribosome_assets import Assetlist, RibosomeAssets, obtain_assets, obtain_assets_threadpool
from ribctl.lib.util_taxonomy import filter_by_parent_tax

# from api.db.ribosomexyz import ribosomexyzDB # from api.rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, RIBETL_DATA
arg = argparse.ArgumentParser(description='RibCtl - A tool to control the ribosome database')

arg.add_argument('-getall', '--obtain_all_structures', action='store_true')
arg.add_argument('-syncall', '--sync_all_structures_with_pdb',   action='store_true')
arg.add_argument('-o', '--obtain', type=str)
arg.add_argument('-ls', '--list_structs', action='store_true')

struct_filter_arggroup = arg.add_argument_group("Structure filtering options")
struct_filter_arggroup.add_argument('-tax', '--taxid', type=int)

arg.add_argument('-ttt', '--test', action='store_true')

args = arg.parse_args()

if args.obtain_all_structures:
    ASL = Assetlist(profile=True)
    obtain_assets_threadpool(
        [],
        ASL,
        workers=16,
        get_all=True,
        # overwrite=True
    )

if args.obtain:
    RCSB_ID = str(args.obtain)
    loop    = asyncio.get_event_loop()
    loop.run_until_complete(
        obtain_assets(
            RCSB_ID,
            Assetlist(profile=True),
            overwrite=True
        )
    )
if args.test:
    i = 0
    nascent_chains = []
    for struct in os.listdir(RIBETL_DATA):
        if struct == '6OXI':
            continue
        if len(struct) != 4:
            continue
        else:
            R = RibosomeAssets(struct)
            s = R.profile()
            for chain in [*(s.polymeric_factors if s.polymeric_factors != None else []), *s.proteins]:
                if 'Nascent Chain' in chain.nomenclature:
                    nascent_chains.append(chain.dict())
        with open('nascent_chains.json', 'w') as f:
            json.dump(nascent_chains, f, indent=4)

if args.list_structs:
    all_structs = os.listdir(RIBETL_DATA)
    if args.taxid:
        print(filter_by_parent_tax(args.taxid))
    else:
        for struct in all_structs:
            print(struct)
            # all_structs = [struct for struct in all_structs if struct.startswith(str(taxid))]



# if args.db:
#     db = ribosomexyzDB(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
#     fire.Fire(db)

# if args.test:
#     #TODO: Make a test suite
#     qo = ribosomexyzDB(NEO4J_URI, NEO4J_USER, NEO4J_PASSWORD)
#     q  = qo.get_all_structures()

#     ss  = qo.get_struct("3J7Z")
#     NeoStruct.validate(ss)

#     ll  = qo.get_individual_ligand("ERY")
#     Ligand.validate(ll)

#     lli = qo.get_all_ligandlike()
#     for lli_i in lli:
#         LigandlikeInstance.validate(lli_i)

#     st = qo.get_RibosomeStructure('3J7Z')
#     RibosomeStructure.validate(st)


#     lig_by_struct = qo.get_ligands_by_struct()
#     # for lbs in lig_by_struct:
#         # LigandsByStruct.validate(lbs)
#         # print("-----------")
#         # pprint(lig_by_struct)

#     spw = qo.match_structs_w_proteins( [ 'uL22','uL4'] )
#     # for struct in spw:
#     #     print(struct)

#     fs  = qo.get_full_structure("5AFI")
#     NeoStruct.validate(fs)

#     bc = qo.get_banclass_for_chain('5AFI', 'I')
#     # pprint(bc)
#     meta = qo.get_banclasses_metadata('e','SSU')

#     nomclassses = qo.list_nom_classes()

#     members = qo.gmo_nom_class('uL2')
#     for m in members:
#         NomenclatureClassMember.validate(m)
#         # print(m)

#     protnum = qo.proteins_number()
#     print(protnum)

#     rbs = qo.get_rnas_by_struct()
#     for r in rbs:
#         ExogenousRNAByStruct.validate(r)

#     rc = qo.get_rna_class('5SrRNA')

#     for rrr in rc:
#         NomenclatureClassMember.validate(rrr)

# url = 'https://rest.uniprot.org/uniprotkb?query=annotation:(type:positional)ANDreviewed:yes&format=fasta&limit=100'


# with requests.get(url.encode()) as request:
#     print(request)
#     request.raise_for_status()
#     with open('rbfa.fasta.gz', 'wb') as f:
#         for chunk in request.iter_content(chunk_size=2**20):
#             print("processing chunk", chunk)
#             f.write(chunk)


# (uniref_cluster_90:UniRef90_Q6F7I0) NOT (accession:A7ZS64
#  curl https://rest.uniprot.org/uniprotkb/stream\?fields\=accession%2Corganism_id%2Csequence\&format\=tsv\&query\=%28%28uniref_cluster_90%3AUniRef90_Q6F7I0%29+NOT+%28accession%3AA7ZS64%29%29         


# import requests
# url = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=(rbfa)AND(reviewed:true)'
# all_fastas = requests.get(url).text
# print(all_fastas)