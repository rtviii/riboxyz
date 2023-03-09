#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
import json
from pprint import pprint

from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.data import QueryOps
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import Ligand, RibosomeStructure
from api.schema.data_requests import LigandsByStruct
from api.schema.v0 import LigandInstance, LigandlikeInstance, NeoStruct
from ribctl.lib.types.types_ribosome_assets import  RibosomeAssets
from ribctl.lib.struct_rcsb_api import current_rcsb_structs

import logging

arg = argparse.ArgumentParser(description='RibCtl - A simple tool to control the ribosome')

arg.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')
arg.add_argument('-ts','--ts_typegen', type=str,nargs='+', help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
arg.add_argument('-s','--structure', type=str)
arg.add_argument('-db','--database', action='store_true')
arg.add_argument('-o','--obtain', type=str)
arg.add_argument('-q','--query', action='store_true')
arg.add_argument('-pdbsync','--sync_rcsb', action='store_true')

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






        



        # try:
        #     # NeoStruct.validate(q)

        # except Exception as e:
        #     print("Failed to validate NeoStruct:", q.struct.rcsb_id)

    # print(qo.driver)
    
    # D        = Neo4jDB()
    # with qo.driver.session() as session:
    #     def _(tx: Transaction | ManagedTransaction):
    #         return tx.run("""//
    #         match (n:RibosomeStructure) return n.rcsb_id
    #         """).values()
    #     print(session.read_transaction(_))


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


