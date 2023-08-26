import argparse

parser = argparse.ArgumentParser(description="Your CLI Description")

# Global argument (affecting all commands)
parser.add_argument('--global-option', help="Global option description")

subparsers = parser.add_subparsers(title='Subcommands', dest='command')

# Command 1
command1_parser = subparsers.add_parser('command1', help='Command 1 help')
command1_parser.add_argument('arg1', help='Argument for command 1')

# Command 2
command2_parser = subparsers.add_parser('command2', help='Command 2 help')
command2_parser.add_argument('--option2', help='Option for command 2')

# Subcommand of command 2
subcommand_parser = command2_parser.add_subparsers(title='Subcommands', dest='subcommand2')
subcommand_parser.add_parser('subcommand2a', help='Subcommand 2a help')
subcommand_parser.add_parser('subcommand2b', help='Subcommand 2b help')

def execute_command1(args):
    print("Executing command 1 with arg:", args.arg1)

def execute_command2(args):
    print("Executing command 2")
    if args.subcommand2 == 'subcommand2a':
        print("Executing subcommand 2a")
    elif args.subcommand2 == 'subcommand2b':
        print("Executing subcommand 2b")
command1_parser.set_defaults(func=execute_command1)
command2_parser.set_defaults(func=execute_command2)



args = parser.parse_args()
if hasattr(args, 'func'):
    args.func(args)
else:
    parser.print_help()
