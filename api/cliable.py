#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
from pprint import pprint
from pydantic2ts import generate_typescript_defs
from api.ribctl.db.driver import Neo4jDB, init_driver
from api.ribctl.lib.types.types_ribosome import RibosomeStructure
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from ribctl.lib.struct_rcsb_api import current_rcsb_structs

arg = argparse.ArgumentParser(description='RibCtl - A simple tool to control the ribosome')

arg.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')
arg.add_argument('-ts','--ts_typegen', type=str,nargs='+', help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
arg.add_argument('-s','--structure', type=str)
arg.add_argument('-db','--database', action='store_true')
arg.add_argument('-o','--obtain', type=str)
arg.add_argument('-pdbsync','--sync_rcsb', action='store_true')

args = arg.parse_args()





if args.sync_rcsb:
    D = Neo4jDB()

    rcsb_structs  = current_rcsb_structs()
    neo4j_structs = D.get_all_structs()
    print(neo4j_structs)



    





if args.obtain:
    rcsb_id = args.obtain
    print(rcsb_id)
    print(rcsb_id)

if args.structure:
   
    r  = RibosomeAssets(args.structure)
    d  = r.json_profile()
    rs = RibosomeStructure(**d)

    pprint(rs)

if args.database:
    D = Neo4jDB()

    for i in ["4UG0","3J9M","5AFI","7K00"]:
        ra = RibosomeAssets(i)
        ra._verify_json_profile(True)
        jp = ra.json_profile()
        pprint(jp)
        # print(D.add_structure(ra))