import sys
import warnings
from Bio import BiopythonExperimentalWarning, BiopythonWarning, BiopythonDeprecationWarning

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonExperimentalWarning)

with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonDeprecationWarning)

from ribctl.lib.schema.types_ribosome import RibosomeStructure
sys.dont_write_bytecode = True
import argparse
import json
from ribctl.cli.ls import cmd_ls
from ribctl.cli.sync import cmd_db
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand
from ribctl.lib.schema.types_binding_site import BindingSite



parser = argparse.ArgumentParser(description="Command line interface for the `ribctl` package.")


parser.add_argument('--verify_schema', action='store_true', help="Verify the schema for every file in the database")

subparsers     = parser.add_subparsers(title='Subcommands', dest='command')
parser_cmd_etl = subparsers.add_parser('etl', help='Acquisition and processing of ribosomal structures and assets.')

parser_lig = subparsers.add_parser('lig', help='ligands')
parser_lig.add_argument('--chemid', type=str, required=True, help='Chemical identifier')
parser_lig.add_argument('--src', type=str, required=True, help='Source file or path')
parser_lig.add_argument('--dest', type=str, required=True, help='Destination file or path')

parser_cmd_etl.add_argument('-getall'      , '--obtain_all_structures', action='store_true')
parser_cmd_etl.add_argument('--rcsb_id'               , dest   ='rcsb_id'    )

parser_cmd_etl.add_argument('-transpose_ligand', dest   ='transpose_ligand'    )
parser_cmd_etl.add_argument('--profile'                 , action ='store_true' )
parser_cmd_etl.add_argument('--ptc_coords'              , action ='store_true' )
parser_cmd_etl.add_argument('--cif'                     , action ='store_true' )
parser_cmd_etl.add_argument('--cif_modified_and_chains' , action ='store_true' )
parser_cmd_etl.add_argument('--ligands'                 , action ='store_true' )
parser_cmd_etl.add_argument('--png_thumbnail'           , action ='store_true' )
parser_cmd_etl.add_argument('--overwrite'               , action ='store_true' )


parser_cmd_etl.add_argument('--ncbi_init' , action ='store_true' )

import asyncio
import os
from ribctl import ASSETS, RIBETL_DATA
from ribctl.etl.etl_obtain import asset_routines, obtain_asssets_threadpool
from ribctl.etl.etl_ribosome_ops import Assetlist, Structure

def cmd_etl(args):

    ASL = Assetlist(
        profile                 = False,
        ptc_coords              = False,
        cif                     = False,
        chains = False,
        ligands                 = False,
        png_thumbnail           = False,
    )

    if args.profile:
        ASL.profile=True

    if args.ptc_coords:
        ASL.ptc_coords=True

    if args.cif:
        ASL.cif=True

    if args.cif_modified_and_chains:
        ASL.chains=True

    if args.ligands:
        ASL.ligands=True

    if args.png_thumbnail:
        ASL.png_thumbnail=True

    #All structures
    if args.obtain_all_structures:
        obtain_asssets_threadpool(
            ASL,
            workers   = 4,
            overwrite = args.overwrite or False
        )

    if args.rcsb_id:
        RCSB_ID = str(args.rcsb_id)
        loop    = asyncio.get_event_loop()
        loop.run_until_complete(
            asset_routines(
                RCSB_ID,
                ASL,
                args.overwrite or False
            )
        )
    else:
        parser_cmd_etl.print_help()


parser_cmd_etl.set_defaults(func=cmd_etl)
# parser_cmd_etl.print_help()

#! -------------------------- sync     -------------------------- #
parser_sync = subparsers.add_parser('db', help='Syncronization with the PDB, updates and database uploads')
parser_sync.set_defaults(func=cmd_db)




# #! -------------------------- ls       -------------------------- #
parser_cmd_ls = subparsers.add_parser('ls', help='List information')
parser_cmd_ls.add_argument('-struct', help="Structure ID")
parser_cmd_ls.add_argument('-taxid', '--taxid', help="Species ID")
parser_cmd_ls.add_argument('-elem', '--subelement', help="Subelement type (rna,protein,ligand)")
parser_cmd_ls.set_defaults(func=cmd_ls)

# #! -------------------------- lig       -------------------------- #
def cmd_lig(args):
    chemid = args.chemid.upper()
    src    = args.src.upper()
    dest   = args.dest.upper()

    
    # with open(BindingSite.path_nonpoly_ligand(chemid, src), 'r') as infile:
    # ligdict      = json.load(infile)
    bsite_src    = BindingSite.parse_file(BindingSite.path_nonpoly_ligand(src,chemid))
    dest_profile = Structure(dest).profile()
    # BindingSite.parse_obj()
    transpose = init_transpose_ligand(dest_profile, bsite_src)
    t         = transpose
    s= {}
    for key in t.dict().keys():
        s[key.value] = t[key]

    dest_path = os.path.join(RIBETL_DATA, dest, f'PREDICTION_{chemid}_{src}_{dest}.json'.upper())
    with open(dest_path, 'w+') as outfile:
        json.dump(s, outfile, indent=4)
        print("saved ", dest_path)
    

parser_lig.set_defaults(func=cmd_lig)

#! -------------------------- Filerts and options -------------------------- #
parser.add_argument('--has_protein', type=parse_comma_separated_list, help="Global option description")
parser.add_argument('--taxid')
parser.add_argument('--verify', action='store_true')
parser.add_argument('--t', action='store_true')
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------


# def verify_structure_profile_schema(rcsb_id:str):
#     s = RibosomeAssets(rcsb_id).profile().model_dump_json()
#     try:
#         RibosomeStructure.model_validate_json(s)
#         return True
#     except Exception as e:
#         print(e)
#         return False

# def verify_profile_exists(rcsb_id:str):
#     return os.path.exists(RibosomeAssets(rcsb_id)._json_profile_filepath())

try:
    args = parser.parse_args()

    if args.verify:
        for struct in os.listdir(RIBETL_DATA):
            print(struct, verify_profile_exists(struct))
            if not verify_profile_exists(struct):
                asyncio.run(asset_routines(struct, Assetlist(profile=True), overwrite=True))
        exit(0)

    # if args.ncbi_init:

    #     for struct in os.listdir(RIBETL_DATA):
    #         print(struct, verify_profile_exists(struct))
    #         if not verify_profile_exists(struct):
    #             asyncio.run(obtain_assets(struct, Assetlist(profile=True), overwrite=True))
    #     exit(0)

    if hasattr(args, 'func'):
        args.func(args)

    elif args.verify_schema:
        all_structs = os.listdir(RIBETL_DATA)
        tally       = { "valid": [], "invalid": []}

        for struct in all_structs:
            if  not os.path.exists(Structure(struct)._json_profile_filepath()):
                continue
            if not verify_structure_profile_schema(struct):
                tally["invalid"].append(struct)
            else:
                tally["valid"].append(struct)
        print("Valid:" , len(tally["valid"]))
        print("Invalid:", len(tally["invalid"]))
    else:
        parser.print_help()
except (BrokenPipeError, IOError) as e:
    print("BrokenPipeError or IOError", e)
    pass




# Notes
# `awk '/ERROR/ {print $3}' | sed 's/:.*$//'`  to get every structure in the log file that failed
