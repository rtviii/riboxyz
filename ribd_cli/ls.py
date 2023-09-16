
import os
from ribctl import RIBETL_DATA


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