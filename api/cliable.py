#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
from ribctl.db.data import QueryOps
from ribctl.db.ribosomexyz import Neo4jDB
from ribctl.lib.types.types_ribosome import ExogenousRNAByStruct, Ligand, ProteinClass, RibosomeStructure
from schema.data_requests import BanclassMetadata, LigandsByStruct
from schema.v0 import BanClassMetadata, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from ribctl.lib.types.types_ribosome_assets import  RibosomeAssets
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
from processing import mp

import logging

arg = argparse.ArgumentParser(description='RibCtl - A simple tool to control the ribosome database')

arg.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')
arg.add_argument('-ts','--ts_typegen', type=str,nargs='+', help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
arg.add_argument('-s','--structure', type=str)
arg.add_argument('-db','--database', action='store_true')
arg.add_argument('-o','--obtain', type=str)
arg.add_argument('-q','--query', action='store_true')
arg.add_argument('-pdbsync','--sync_rcsb', action='store_true')
arg.add_argument('-p','--process', action='store_true')

args = arg.parse_args()

logging.basicConfig(level=logging.ERROR, filemode='w', format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

if args.sync_rcsb:

    D             = Neo4jDB()
    rcsb_structs  = current_rcsb_structs()
    neo4j_structs = D.get_all_structs()
    print(neo4j_structs)

if args.obtain:
    rcsb_id = args.obtain

if args.structure:
    r  = RibosomeAssets(args.structure)
    d  = r.json_profile()
    rs = RibosomeStructure(**d)

if args.query:

    #TODO: Make a test suite
    qo = QueryOps()
    q  = qo.get_all_structures()

    ss  = qo.get_struct("3J7Z")
    NeoStruct.validate(ss)

    ll  = qo.get_individual_ligand("ERY")
    Ligand.validate(ll)

    lli = qo.get_all_ligandlike()
    for lli_i in lli:
        LigandlikeInstance.validate(lli_i)
        # pprint(lli_i)


    st = qo.get_RibosomeStructure('3J7Z')
    RibosomeStructure.validate(st)


    lig_by_struct = qo.get_ligands_by_struct()
    # for lbs in lig_by_struct:
        # LigandsByStruct.validate(lbs)
        # print("-----------")
        # pprint(lig_by_struct)

    spw = qo.match_structs_w_proteins( [ 'uL22','uL4'] )
    # for struct in spw:
    #     print(struct)

    fs  = qo.get_full_structure("5AFI")
    NeoStruct.validate(fs)

    bc = qo.get_banclass_for_chain('5AFI', 'I')
    # pprint(bc)
    meta = qo.get_banclasses_metadata('e','SSU')

    nomclassses = qo.list_nom_classes()

    members = qo.gmo_nom_class('uL2')
    for m in members:
        NomenclatureClassMember.validate(m)
        # print(m)

    protnum = qo.proteins_number()
    print(protnum)

    rbs = qo.get_rnas_by_struct()
    for r in rbs:
        ExogenousRNAByStruct.validate(r)

    rc = qo.get_rna_class('5SrRNA')

    for rrr in rc:
        NomenclatureClassMember.validate(rrr)

if args.database:

    D        = Neo4jDB()
    synced   = D.get_all_structs()
    unsynced = sorted(current_rcsb_structs())

    for rcsb_id in unsynced:
        assets = RibosomeAssets(rcsb_id)
        try:
            assets._verify_json_profile(True)
            D.add_structure(assets)

        except Exception as e:
            print(e)
            logger.error("Exception occurred:", exc_info=True)

if args.process:
    mp()