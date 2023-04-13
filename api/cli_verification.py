# TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.
# export DJANGO_SETTINGS_MODULE=api.rbxz_bend                                                                                         [cli]
# import argparse
from pprint import pprint
from api.rbxz_bend.db.ribosomexyz import ribosomexyzDB
from api.rbxz_bend.settings import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from ribctl.lib.types.types_ribosome import Ligand, ProteinClass, RibosomeStructure
from schema.data_requests import BanclassMetadata, LigandsByStruct
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from api.ribctl.lib.types.ribosome_assets import RibosomeAssets
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
import fire

# import logging

# arg = argparse.ArgumentParser(
#     description='RibCtl - A tool to control the ribosome database')

# arg.add_argument('-dbname', '--database_name', type=str)
# arg.add_argument('-s', '--structure', type=str)
# arg.add_argument('-o', '--obtain', type=str)
# arg.add_argument('-ttt', '--test', action='store_true')
# # arg.add_argument('-pdbsync', '--sync_rcsb' , action='store_true'                                                                                                         )
# args = arg.parse_args()

# logging.basicConfig(level=logging.ERROR, filemode='w',
#                     format='%(asctime)s - %(levelname)s - %(message)s')
# logger = logging.getLogger(__name__)


# if args.obtain:
#     rcsb_id = args.obtain
# db.get_full_structure('3J7Z')
# if args.structure:
#     r = RibosomeAssets(args.structure)
#     d = r.json_profile()
#     rs = RibosomeStructure.parse_obj(d)

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
