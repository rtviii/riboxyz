import argparse
import json
from pprint import pprint

# from ribctl.cli.etl import cmd_etl
from ribctl.cli.ls import cmd_ls
from ribctl.cli.sync import cmd_sync
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand
from ribctl.lib.ribosome_types.types_binding_site import BindingSite

def parse_comma_separated_list(value):
    return value.split(',')






parser     = argparse.ArgumentParser(description="Command line interface for the `ribctl` package.")
subparsers = parser.add_subparsers(title='Subcommands', dest='command')

#! -------------------------- --- -------------------------- #
#! -------------------------- fasta -------------------------- #
#! -------------------------- --- -------------------------- #

parser_cmd_fasta = subparsers.add_parser('fasta', help='Dealing with sequences')

#! -------------------------- --- -------------------------- #
#! -------------------------- etl -------------------------- #
#! -------------------------- --- -------------------------- #

parser_cmd_etl = subparsers.add_parser('etl', help='Acquisition and processing of ribosomal structures and assets.')


parser_lig = subparsers.add_parser('lig', help='ligands')
parser_lig.add_argument('--chemid', type=str, required=True, help='Chemical identifier')
parser_lig.add_argument('--src', type=str, required=True, help='Source file or path')
parser_lig.add_argument('--dest', type=str, required=True, help='Destination file or path')




parser_cmd_etl.add_argument('-getall'      , '--obtain_all_structures', action='store_true')
parser_cmd_etl.add_argument('-struct'               , dest   ='rcsb_id'    )

parser_cmd_etl.add_argument('-transpose_ligand', dest   ='transpose_ligand'    )
parser_cmd_etl.add_argument('--profile'                 , action ='store_true' )
parser_cmd_etl.add_argument('--ptc_coords'              , action ='store_true' )
parser_cmd_etl.add_argument('--cif'                     , action ='store_true' )
parser_cmd_etl.add_argument('--cif_modified_and_chains' , action ='store_true' )
parser_cmd_etl.add_argument('--ligands'                 , action ='store_true' )
parser_cmd_etl.add_argument('--png_thumbnail'           , action ='store_true' )
parser_cmd_etl.add_argument('--overwrite'               , action ='store_true' )



import asyncio
import os
from ribctl import ASSETS, RIBETL_DATA
from ribctl.etl.obtain import obtain_assets, obtain_assets_threadpool
from ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets

def cmd_etl(args):

    ASL = Assetlist()

    if args.profile:
        ASL.profile=True

    if args.ptc_coords:
        ASL.ptc_coords=True

    if args.cif:
        ASL.cif=True

    if args.cif_modified_and_chains:
        ASL.cif_modified_and_chains=True

    if args.ligands:
        ASL.ligands=True

    if args.png_thumbnail:
        ASL.png_thumbnail=True

    #All structures
    if args.obtain_all_structures:
        obtain_assets_threadpool(
            [],
            ASL,
            workers=4,
            get_all=True,
            overwrite=args.overwrite or False
        )

    if args.rcsb_id:
        RCSB_ID = str(args.rcsb_id)
        loop    = asyncio.get_event_loop()
        loop.run_until_complete(
            obtain_assets(
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
parser_sync = subparsers.add_parser('sync_db', help='Syncronization with the PDB, updates and database uploads')
# parser_cmd2sub = parser_sync.add_subparsers(title='Subcommands', dest='subcommand2')
# parser_cmd2sub.add_parser('db', help='Upload local structures to the neo4j database')
parser_sync.set_defaults(func=cmd_sync)



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
    dest_profile = RibosomeAssets(dest).profile()
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
parser.add_argument('--t', action='store_true')
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------







args = parser.parse_args()
print(args)
if args.t:
    ...
    # test()
else:
    if hasattr(args, 'func'):
        args.func(args)

    else:
        parser.print_help()

# Notes

# `awk '/ERROR/ {print $3}' | sed 's/:.*$//'`  to get every structure in the log file that failed
