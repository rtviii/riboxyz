#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
from pprint import pprint
from pydantic2ts import generate_typescript_defs
from api.ribctl.db.driver import Neo4jDB, init_driver
from api.ribctl.db import constraints, structure, rna, proteins
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets

arg = argparse.ArgumentParser(description='RibCtl - A simple tool to control the ribosome')

arg.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')
arg.add_argument('-ts','--ts_typegen', type=str,nargs='+', help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
arg.add_argument('-s','--structure', type=str)
arg.add_argument('-db','--database', action='store_true')
args = arg.parse_args()


if args.structure:
    RibosomeAssets(args.structure)

if args.database:
    D = Neo4jDB()


    D.see_constraints()
    D.init_constraints()
    pprint(D.see_constraints())

