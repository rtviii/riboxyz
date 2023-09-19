
import json
import os
from pprint import pprint
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets


def cmd_ls(args):

    all_structs = os.listdir(RIBETL_DATA)
    if args.struct != None:

        if "." in args.struct:
            rcsb_id, auth_asym_id = args.struct.split(".")
        rp = RibosomeAssets(rcsb_id)
        if "." in args.struct:
            chain, rp_class = rp.get_chain_by_auth_asym_id(auth_asym_id)
            pprint(json.loads(chain.json()))
        else:
            pprint(json.loads(RibosomeAssets(args.struct).profile().json()))


    elif args.species != None:
        print("Listing species information for", args.species)
    elif args.subelement != None:
        print("Listing subelement information for", args.subelement)
    else:
        print("Listing all information")