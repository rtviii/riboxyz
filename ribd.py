import argparse
from ribd_cli.etl import cmd_etl
from ribd_cli.ls import cmd_ls
from ribd_cli.sync import cmd_sync

def parse_comma_separated_list(value):
    return value.split(',')



parser     = argparse.ArgumentParser(description="Command line interface for the `ribctl` package.")
subparsers = parser.add_subparsers(title='Subcommands', dest='command')

#! -------------------------- --- -------------------------- #
#! -------------------------- etl -------------------------- #
#! -------------------------- --- -------------------------- #

parser_cmd_etl = subparsers.add_parser('etl', help='Acquisition and processing of ribosomal structures and assets.')

parser_cmd_etl.add_argument('-getall'      , '--obtain_all_structures', action='store_true')
parser_cmd_etl.add_argument('-struct'               , dest   ='rcsb_id'    )

parser_cmd_etl.add_argument('--profile'                 , action ='store_true' )
parser_cmd_etl.add_argument('--ptc_coords'              , action ='store_true' )
parser_cmd_etl.add_argument('--cif'                     , action ='store_true' )
parser_cmd_etl.add_argument('--cif_modified_and_chains' , action ='store_true' )
parser_cmd_etl.add_argument('--factors_and_ligands'     , action ='store_true' )
parser_cmd_etl.add_argument('--png_thumbnail'           , action ='store_true' )
parser_cmd_etl.add_argument('--overwrite'               , action ='store_true' )
parser_cmd_etl.set_defaults(func=cmd_etl)


#! -------------------------- -------- -------------------------- #
#! -------------------------- sync     -------------------------- #
#! -------------------------- -------- -------------------------- #
parser_sync = subparsers.add_parser('sync_db', help='Syncronization with the PDB, updates and database uploads')
# parser_cmd2sub = parser_sync.add_subparsers(title='Subcommands', dest='subcommand2')
# parser_cmd2sub.add_parser('db', help='Upload local structures to the neo4j database')
parser_sync.set_defaults(func=cmd_sync)


# #! -------------------------- -------- -------------------------- #
# #! -------------------------- ls       -------------------------- #
# #! -------------------------- -------- -------------------------- #
parser_cmd_ls = subparsers.add_parser('ls', help='List information')
parser_cmd_ls.add_argument('-struct', help="Structure ID")
parser_cmd_ls.add_argument('-spec', '--species', help="Species ID")
parser_cmd_ls.add_argument('-elem', '--subelement', help="Subelement type (rna,protein,ligand)")
parser_cmd_ls.set_defaults(func=cmd_ls)



#! -------------------------- Filerts and options -------------------------- #
parser.add_argument('--has_protein', type=parse_comma_separated_list, help="Global option description")
parser.add_argument('--taxid')
parser.add_argument('--t', action='store_true')
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------
#?---------------------------------------------------------------------------------------------------------







args = parser.parse_args()
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
