# TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.
# export DJANGO_SETTINGS_MODULE=api.rbxz_bend                                                                                         [cli]
# import argparse
import argparse
from api.ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets, obtain_assets, obtain_assets_processpool, obtain_assets_threadpool 
from api.ribctl.lib.types.types_ribosome import RibosomeStructure
from logs.loggers import updates_logger
from ribctl.etl.struct_rcsb_api import gql_monolith,query_rcsb_api, process_pdb_record

arg = argparse.ArgumentParser(description='RibCtl - A tool to control the ribosome database')

arg.add_argument('-dbname', '--database_name', type=str)
arg.add_argument('-s', '--structure', type=str)
arg.add_argument('-o', '--obtain', type=str)
arg.add_argument('-ttt', '--test', action='store_true')

args = arg.parse_args()

if args.obtain:

    rcsb_list = ['3J7Z', '5afi','4ug0']
    obtain_assets_threadpool(
                        rcsb_list,
                        Assetlist(profile=True,  factors_and_ligands=True),
                        workers=16,
                        get_all=True,
                        replace=True
                      )
if args.structure:
    str = args.structure.upper()
    try:

        struct = process_pdb_record(str)
        RibosomeStructure.parse_obj(struct)
        assets = RibosomeAssets(str)
        assets.save_json_profile(assets._json_profile_filepath(), struct.dict())
        updates_logger.info(f"Saved {str}.json to {assets._json_profile_filepath()}")
    except Exception as e:
        updates_logger.exception(e)

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
