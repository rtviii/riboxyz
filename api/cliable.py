#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
from pprint import pprint

from neo4j import ManagedTransaction, Transaction
from api.ribctl.db.data import QueryOps
from api.ribctl.db.inits.driver import Neo4jDB
from api.ribctl.lib.types.types_ribosome import RibosomeStructure
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
    print(qo.get_all_structures())
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


