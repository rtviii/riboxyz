
import json
import os
from pprint import pprint
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import RibosomeStructure
from ribctl.lib.util_taxonomy import  get_descendants_of, taxid_is_descendant_of,descendants_of_taxid


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

    elif args.taxid != None:
        print("Listing species information for", args.taxid)
        all_structs = os.listdir(RIBETL_DATA)
        pdbid_taxid_tuples:list = []    

        for struct in all_structs:
            rp = RibosomeAssets(struct).profile()
            pdbid_taxid_tuples.append(( rp.rcsb_id, rp.src_organism_ids[0] ))

        print(descendants_of_taxid( pdbid_taxid_tuples, int(args.taxid)))

    elif args.subelement != None:
        print("Listing subelement information for", args.subelement)
    else:
        print(all_structs)