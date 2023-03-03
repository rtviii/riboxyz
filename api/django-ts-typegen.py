# TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.
from pydantic2ts import generate_typescript_defs
import argparse
import os
os.environ["DJANGO_SETTINGS_MODULE"] = 'rbxz_bend.settings'

arg = argparse.ArgumentParser(
    description='RibCtl - A simple tool to control the ribosome')
arg.add_argument('-ts', '--ts_typegen', type=str, nargs='+',
                 help='Generate TypeScript types for the API. Point at "types" directory and specify the out directory')
args = arg.parse_args()

DEFAULT_DJANGO_TYPES_DIR = os.path.join(os.path.dirname(__file__), 'schema')
DEFAULT_DJANGO_OUT_DIR = os.path.join(
    os.path.dirname(__file__), '__typescript_d_ts')


def typegen(src, dest):
    for file_in in os.listdir(src):
        if file_in.endswith(".d.py"):
            file_in = os.path.join(src, file_in)
            file_out = os.path.join(dest, file_in.replace(".d.py", ".d.ts"))
            generate_typescript_defs(file_in, file_out)
            print("Generated: {}".format(file_out))


if __name__ == "__main__":
    typegen(DEFAULT_DJANGO_TYPES_DIR, DEFAULT_DJANGO_OUT_DIR)
    exit(0)

if args.ts_typegen:
    if len(args.ts_typegen) != 2:
        raise Exception(
            "You must specify the input and output directories for the TypeScript type generation.")
    [src, dest] = args.ts_typegen
    typegen(src, dest)
