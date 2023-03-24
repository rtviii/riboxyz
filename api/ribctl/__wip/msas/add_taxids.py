import os
from pprint import pprint
import subprocess
import sys
import numpy as np
import prody as prd
from prody import MSA, Sequence
import requests
prd.confProDy(verbosity='none')

def muscle_combine_profiles(msa_path1: str, msa_path2: str, out_filepath: str):
    """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
    cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]

    subprocess.Popen(cmd,
                      stdout=subprocess.PIPE,
                      stderr=subprocess.PIPE,
                      env=os.environ.copy()
                      ).wait()

    sys.stdout.flush()

def barr2str (bArr):
    return ''.join([ x.decode("utf-8") for x in bArr])
    
infer_subunit = lambda polymer_class : 'LSU' if 'L' in polymer_class else 'SSU'
datapath      = lambda subunit, polymer_class : '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/{}/{}_ribovision.fasta'.format(subunit, polymer_class)


nomclass      = 'uL23'
path          = datapath(infer_subunit(nomclass),nomclass)



def msa_add_taxonomic_ids(msa_path:str):
    msa_main:MSA = prd.parseMSA(msa_path)
    url = lambda protein_id: f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id={protein_id}&retmode=json"
    seq:Sequence

    _sequences_new = []
    _labels_new   = []

    for seq in msa_main:
        original_label  = seq.getLabel()
        nih_protein_id  = original_label.split('|')[1]
        sequence_proper = barr2str(seq.getArray())

        response        = requests.get(url(nih_protein_id))
        if response.status_code == 200:
            protein_summary = response.json()
            taxid           = protein_summary["result"][nih_protein_id]["taxid"]
            new_label = f"{original_label}|{taxid}"

            _sequences_new.append(sequence_proper)
            _labels_new.append(new_label)
        else:
            print("Omtitting {} Error: ".format(original_label), response.status_code)

    prd.writeMSA(msa_path, MSA(np.array(_sequences_new),nomclass,labels=_labels_new))






# star -----------------------------
# from ete3 import NcbiTaxa

# # create an instance of the NcbiTaxa class
# ncbi = NcbiTaxa()

# # define the tax ID of interest
# tax_id = 9606  # for human

# # retrieve the lineage
# lineage = ncbi.get_lineage(tax_id)

# # retrieve the scientific names of the lineage
# lineage_names = ncbi.get_taxid_translator(lineage)

# # print the lineage
# for taxid, name in lineage_names.items():
#     print(f"{taxid}\t{name}")
