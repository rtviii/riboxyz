# TODO: ( Bring up ) Lift all the useful top-level functions from the package to this level.
import os
os.environ["DJANGO_SETTINGS_MODULE"] = 'rbxz_bend.settings'
from pydantic2ts import generate_typescript_defs
import argparse


DEFAULT_DJANGO_TYPES_DIR = os.path.join(os.path.dirname(__file__), 'schema')
DEFAULT_DJANGO_OUT_DIR   = os.path.join(os.path.dirname(__file__), '__typescript_d_ts')

def typegen(schema_dir, d_ts_dir):
    if not os.path.exists(d_ts_dir):
        os.mkdir(d_ts_dir)

    for file_in in os.listdir(schema_dir):
        if file_in.endswith(".py"):
            file_in  = os.path.join(schema_dir, file_in)
            # file_out = '{}/{}.d.ts'.format(d_ts_dir,os.path.splitext(os.path.basename(file_in))[0])
            file_out = os.path.join(d_ts_dir,os.path.splitext(os.path.basename(file_in))[0] + ".d.ts")
            generate_typescript_defs(file_in, file_out)
            print("[TYPEGEN] Generated: {}".format(file_out))


arg = argparse.ArgumentParser()
arg.add_argument('-ts', '--ts_typegen', type=str, nargs='+',help=' Point at "types" directory and specify the out directory')
args = arg.parse_args()

if __name__ == "__main__":
    print("Running typegen")
    typegen(DEFAULT_DJANGO_TYPES_DIR, DEFAULT_DJANGO_OUT_DIR)
    exit(0)

if args.ts_typegen:
    if len(args.ts_typegen) != 2:
        raise Exception(
            "You must specify the input and output directories for the TypeScript type generation.")
    [schema_dir, d_ts_dir] = args.ts_typegen
    typegen(schema_dir, d_ts_dir)
