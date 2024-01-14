import json
import os
from pprint import pprint
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libmsa import Fasta
from io import StringIO
import os
from Bio import SeqRecord
import os
import subprocess
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
records = []

for s in RibosomeAssets.list_all_structs():
    try:
        with open(os.path.join(RIBETL_DATA,s,"{}.json".format(s)), 'r') as f:
            struct = json.load(f)
            if 'polymeric_factors' in struct:
                print(True)
                for pf in struct['polymeric_factors']:
                    if 'tRNA' in pf['nomenclature']:
                       seq   = pf['entity_poly_seq_one_letter_code_can']
                       taxid = pf['src_organism_ids'][0]
                       desc  = f"{ pf['parent_rcsb_id'] }.{pf['auth_asym_id']} {pf['rcsb_pdbx_description']}"
                       if "hydrolase" in desc:
                           print("skipping hydrolase")
                           continue
                       if len(seq) < 50:
                           print("skipping aberrant")
                           continue


                       s     = SeqRecord(Seq(seq), id=str(taxid), description=desc)
                       records.append(s)
    except Exception as e:
        print(e)
Fasta.write_fasta(records, "/home/rtviii/dev/riboxyz/ribctl/assets/fasta_trna/trna.fasta")
