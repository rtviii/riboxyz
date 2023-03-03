#TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.

import argparse
import os
from pydantic2ts import generate_typescript_defs

from ribctl.lib.types.types_ribosome_assets import RibosomeAssets

arg = argparse.ArgumentParser(description='RibCtl - A simple tool to control the ribosome')

arg.add_argument('-v', '--version', action='version', version='%(prog)s 0.1.0')
arg.add_argument('-ts','--ts_typegen', type=str,nargs='+', help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
arg.add_argument('-s','--structure', type=str)

args = arg.parse_args()

if args.ts_typegen:
    s = args.ts_typegen
    generate_typescript_defs(s[0], s[1])


if args.structure:
    RibosomeAssets(args.structure)
    print(args.structure)