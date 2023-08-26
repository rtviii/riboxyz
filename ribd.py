import argparse

parser = argparse.ArgumentParser(description="Command line interface for the `ribctl` package.")




# arg.add_argument('-getall', '--obtain_all_structures', action='store_true')
# arg.add_argument('-syncall', '--sync_all_structures_with_pdb',   action='store_true')
# arg.add_argument('-o', '--obtain', type=str)
# arg.add_argument('-ls', '--list_structs', action='store_true')

# struct_filter_arggroup = arg.add_argument_group("Structure filtering options")
# struct_filter_arggroup.add_argument('-tax', '--taxid', type=int)

# arg.add_argument('-ttt', '--test', action='store_true')





subparsers = parser.add_subparsers(title='Subcommands', dest='command')

#! -------------------------- --- -------------------------- #
#! -------------------------- ETL -------------------------- #
#! -------------------------- --- -------------------------- #

parser_cmd_etl = subparsers.add_parser('etl', help='Acquisition and processing of ribosomal structures and assets.')
parser_cmd_etl.add_argument('-getall', '--obtain_all_structures', action='store_true')

def cmd_etl(args):
    print("Executing command 1 with arg:", args.obtain_all_structures)

parser.set_defaults(func=cmd_etl)


#! -------------------------- -------- -------------------------- #
#! -------------------------- Command2 -------------------------- #
#! -------------------------- -------- -------------------------- #
parser_cmd2 = subparsers.add_parser('sync', help='Syncronization with the PDB, updates and database uploads')
parser_cmd2.add_argument('--option2', help='Option for command 2')

parser_cmd2sub = parser_cmd2.add_subparsers(title='Subcommands', dest='subcommand2')

parser_cmd2sub.add_parser('subcommand2a', help='Subcommand 2a help')
parser_cmd2sub.add_parser('subcommand2b', help='Subcommand 2b help')

def cmd_2(args):
    print("Executing command 2")
    if args.subcommand2 == 'subcommand2a':
        print("Executing subcommand 2a")
    elif args.subcommand2 == 'subcommand2b':
        print("Executing subcommand 2b")

parser_cmd2.set_defaults(func=cmd_2)


#! -------------------------- -------- -------------------------- #
#! -------------------------- Command2 -------------------------- #
#! -------------------------- -------- -------------------------- #

parser_cmd2 = subparsers.add_parser('sync', help='Syncronization with the PDB, updates and database uploads')
parser_cmd2.add_argument('--option2', help='Option for command 2')

parser_cmd2sub = parser_cmd2.add_subparsers(title='Subcommands', dest='subcommand2')

parser_cmd2sub.add_parser('subcommand2a', help='Subcommand 2a help')
parser_cmd2sub.add_parser('subcommand2b', help='Subcommand 2b help')

def cmd_2(args):
    print("Executing command 2")
    if args.subcommand2 == 'subcommand2a':
        print("Executing subcommand 2a")
    elif args.subcommand2 == 'subcommand2b':
        print("Executing subcommand 2b")

parser_cmd2.set_defaults(func=cmd_2)

#! -------------------------- ------------------- -------------------------- #
#! -------------------------- Filerts and options -------------------------- #
#! -------------------------- ------------------- -------------------------- #
parser.add_argument('--modifier-option', help="Global option description")


args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
    print(args.modifier_option, "++")
else:
    parser.print_help()
