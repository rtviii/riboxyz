from functools import reduce
from io import StringIO
from itertools import tee
import json
import logging
import os
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Iterator
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.classification import classify_sequence, classify_subchains
from ribctl.lib.msalib import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.lib.ribosome_types.types_ribosome import ProteinClass, ProteinClassEnum, RibosomeStructure
from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
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
        struct, _ = f.split(".")
        # print(struct)
        rib_asset = RibosomeAssets(struct)
        # print(json_repr)


        filepath =rib_asset._json_profile_filepath()
        profile = rib_asset.profile().json()




        # rib_asset.write_own_json_profile()
        # print(rib_assest.json())

        # with open(os.path.join(nomv2dict, f),'r') as infile:
        #     nomcl = json.load(infile)
        #     pprint(nomcl)

        # with open(os.path.join(nomv2dict, f),'r') as infile:
        #     nomcl = json.load(infile)
        #     pprint(nomcl)
        exit()

    
elif sys.argv[1] =="tunnel":
    print("tunnel")