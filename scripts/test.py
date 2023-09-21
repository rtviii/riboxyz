from functools import reduce
from io import StringIO
from itertools import tee
import json
import logging
import os
from pprint import pprint
from tempfile import NamedTemporaryFile
from typing import Iterator
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from ribctl import ASSETS
from ribctl.lib.classification import classify_sequence, classify_subchains
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
hmm_cachedir = ASSETS['__hmm_cache']
import sys

logger       = logging.getLogger(__name__)
file_handler = logging.FileHandler('classification.log')
log_format   = logging.Formatter('%(asctime)s [%(levelname)s] [%(name)s] %(message)s')
file_handler.setFormatter(log_format)
logger.addHandler(file_handler)

nomv2dict= '/home/rtviii/dev/riboxyz/nomv2'

if sys.argv[1]    == "process_struct":
   rcsb_id         = sys.argv[2].upper()
   rib             = RibosomeAssets(rcsb_id).profile()
   organism_taxid  = rib.src_organism_ids[0]
   prots           = rib.proteins
   result          = classify_subchains(prots)
   print(result)
elif sys.argv[1] == "process_all":
    for rcsb_id in RibosomeAssets.list_all_structs():
        logger.debug("Processing {}".format(rcsb_id))
        rib            = RibosomeAssets(rcsb_id).profile()
        organism_taxid = rib.src_organism_ids[0]
        prots          = rib.proteins
        results        = classify_subchains(prots)
        with open('/home/rtviii/dev/riboxyz/nomv2/{}.json'.format(rcsb_id), 'w') as f:
            json.dump(results, f)


elif sys.argv[1] =="merge_nomenclature":


    for f in os.listdir(nomv2dict):
        
        print(f)
        with open(os.path.join(nomv2dict, f),'r') as infile:
            nomv2=json.load(infile)
        struct, _ = f.split(".")
        rib_asset = RibosomeAssets(struct)


        filepath = rib_asset._json_profile_filepath()
        profile  = rib_asset.profile()
        
        for chain in profile.proteins:
            if chain.auth_asym_id in nomv2:
                if nomv2[chain.auth_asym_id] != None and nomv2[chain.auth_asym_id] not in chain.nomenclature:
                    chain.nomenclature = [nomv2[chain.auth_asym_id]]
                elif nomv2[chain.auth_asym_id] == None and chain.nomenclature != []:
                    chain.nomenclature = []

        rib_asset.write_own_json_profile(new_profile=json.loads(profile.json()), overwrite=True)


elif sys.argv[1] =="tunnel":

    def list_euk_structs():
        EUK_STRUCTS = []
        with open("eukarya_03_07_2023.txt", "r") as data_file:
            for line in data_file:
                structs = line.split(",")
                EUK_STRUCTS = [*EUK_STRUCTS, *structs]
        return EUK_STRUCTS

    if __name__ == "__main__":

        EUK        = list_euk_structs()
        PTC_COORDS = {}
        for RCSB_ID in EUK :
            try:
                print("Processing {}".format(RCSB_ID))
                r= RibosomeAssets(RCSB_ID).profile()
                ress, auth_asym_id = ptc_resdiues_get(RCSB_ID, r.rnas)
                midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)

                residue_labels = [(res.get_resname(), res.id[1]) for res in ress]
                print(residue_labels)

                writeout = {
                    "site_9_residues": [
                        (res.get_resname(), res.get_segid()) for res in ress
                    ],
                    "LSU_rRNA_auth_asym_id": auth_asym_id,
                    "midpoint_coordinates": midpoint_coords,
                }

                PTC_COORDS = {**PTC_COORDS, RCSB_ID: writeout}

            except Exception as e:
                print(e)
   
    print("tunnel")