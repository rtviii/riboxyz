import argparse
import asyncio
import os
from driver import test
from ribctl import RIBETL_DATA
from ribctl.ribosome_assets import Assetlist, obtain_assets, obtain_assets_threadpool




parser     = argparse.ArgumentParser(description="Command line interface for the `ribctl` package.")
subparsers = parser.add_subparsers(title='Subcommands', dest='command')

#! -------------------------- --- -------------------------- #
#! -------------------------- etl -------------------------- #
#! -------------------------- --- -------------------------- #

parser_cmd_etl = subparsers.add_parser('etl', help='Acquisition and processing of ribosomal structures and assets.')
parser_cmd_etl.add_argument('-getall', '--obtain_all_structures', action='store_true')
parser_cmd_etl.add_argument('-get', '--obtain_structure',dest='rcsb_id', help="Get a particular structure by its RCSB_ID")

def cmd_etl(args):
    print("")
    if args.obtain_all_structures:
        ASL = Assetlist(profile=True)
        obtain_assets_threadpool(
            [],
            ASL,
            workers=16,
            get_all=True,
            overwrite=True
        )

    if args.obtain_structure:
        RCSB_ID = str(args.obtain_structure)
        loop    = asyncio.get_event_loop()
        loop.run_until_complete(
            obtain_assets(
                RCSB_ID,
                Assetlist(profile=True),
                overwrite=True
            )
        )

parser_cmd_etl.set_defaults(func=cmd_etl)


#! -------------------------- -------- -------------------------- #
#! -------------------------- sync     -------------------------- #
#! -------------------------- -------- -------------------------- #
parser_sync = subparsers.add_parser('sync_db', help='Syncronization with the PDB, updates and database uploads')
# parser_cmd2sub = parser_sync.add_subparsers(title='Subcommands', dest='subcommand2')
# parser_cmd2sub.add_parser('db', help='Upload local structures to the neo4j database')

def cmd_sync(args):
    print("Uploading structures to Neo4j")

parser_sync.set_defaults(func=cmd_sync)


# #! -------------------------- -------- -------------------------- #
# #! -------------------------- ls       -------------------------- #
# #! -------------------------- -------- -------------------------- #
parser_cmd_ls = subparsers.add_parser('ls', help='List information')

parser_cmd_ls.add_argument('-struct', help="Structure ID")
parser_cmd_ls.add_argument('-spec', '--species', help="Species ID")
parser_cmd_ls.add_argument('-elem', '--subelement', help="Subelement type (rna,protein,ligand)")

def cmd_ls(args):

    all_structs = os.listdir(RIBETL_DATA)

    if args.struct != None:
        print("Listing structure information for", args.struct)
    elif args.species != None:
        print("Listing species information for", args.species)
    elif args.subelement != None:
        print("Listing subelement information for", args.subelement)
    else:
        print("Listing all information")

parser_cmd_ls.set_defaults(func=cmd_ls)


def parse_comma_separated_list(value):
    return value.split(',')

#! -------------------------- Filerts and options -------------------------- #
parser.add_argument('--has_protein', type=parse_comma_separated_list, help="Global option description")
parser.add_argument('--taxid')
parser.add_argument('--t', action='store_true')
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------







args = parser.parse_args()
if args.t:
    test()
else:
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()

# Notes

# `awk '/ERROR/ {print $3}' | sed 's/:.*$//'`  to get every structure in the log file that failed
