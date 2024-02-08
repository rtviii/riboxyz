
import json
import os
from pprint import pprint
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import PolynucleotideClass, PolypeptideClass, RibosomeStructure
# from ribctl.lib.libmsa import  Taxid


def cmd_ls(args):
    all_structs = os.listdir(RIBETL_DATA)

    if args.struct != None:
        rcsb_id, auth_asym_id = None,None
        if "." in args.struct:
            rcsb_id, auth_asym_id = args.struct.split(".")
        else:
            rcsb_id = args.struct


        ribosome_Assets = RibosomeAssets(rcsb_id)

        if "." in args.struct:
            chain, rp_class = ribosome_Assets.get_chain_by_auth_asym_id(auth_asym_id)
            if chain != None:
                print(json.loads(chain.model_dump_json()))
        else:
            print(RibosomeAssets(args.struct).profile().model_dump_json())

    elif args.taxid != None:
        print("Listing species information for", args.taxid)
        all_structs = os.listdir(RIBETL_DATA)
        pdbid_taxid_tuples:list = []    

        for struct in all_structs:
            ribosome_Assets = RibosomeAssets(struct).profile()
            pdbid_taxid_tuples.append(( ribosome_Assets.rcsb_id, ribosome_Assets.src_organism_ids[0] ))

        pprint(Taxid.descendants_of_taxid( pdbid_taxid_tuples, int(args.taxid)))

    elif args.subelement != None:
        subelem = args.subelement
        found   =  []

        try:
            assert(subelem in [_.value for _ in [*list(PolynucleotideClass), *list(PolypeptideClass)]])
        except AssertionError:
            print("Subelement must be one of the following:")
            print([_.value for _ in [*list(PolynucleotideClass), *list(PolypeptideClass)]])
            exit(1)
        for struct in all_structs:
            ra   = RibosomeAssets(struct)
            elem = ra.get_chain_by_polymer_class(subelem)

            if elem  != None:
                found.append(elem)
            else:
                ...

        with open('found_{}.json'.format(subelem), 'w') as outfile:
            json.dump([json.loads(_.json()) for _ in found], outfile, indent=4)
            print("Saved:", 'found_{}.json'.format(subelem))
                
               

           



    else:
        print(all_structs)