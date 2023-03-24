import argparse
import os
from pprint import pprint
import subprocess
import sys
import numpy as np
import prody as prd
from prody import MSA, Sequence, calcShannonEntropy
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
    
# muscle_combine_profiles(chief_msa_path, control_path, '{}with_{}.fasta'.format(chief_msa_path,control))
# h, hwctl = compare_entropies(chief_msa_path, '{}with_{}.fasta'.format(chief_msa_path,control))

# ※ ---------------------
# Goals are to help Shiqi find the the common landmark for the exit port:
    # For every domain:
# 1. Pick the target in one structure X (let's say it's (nucleotide U61 in chain ??) or (amino-acid ?? in uL23)). Call this the origin site: chain + nucleotide/amino-acid.
# 2. Given the species of X, get an alignment of all the other species in this domain. 
# 3. Align the original chain into this MSA (profile-profile in muscle). Call the resulting [columns] aligned to the origin site in the MSA: landmark columns/indices
        # For all other structures of interest in this domain:
        # 4. Get the chain of the same class and align it into the MSA.
        # 5. Track back the reference columns from the MSA-stretched chain to the original (count and exclude the gaps basically)



# ! All the alignments in the proteovision folder are taken for as many domains(2,2157,2759) as possible (at most -- whole 3 domains).
nomclass  = 'uL23'
datapath  = lambda subunit,polymer_class : '/home/rxz/dev/docker_ribxz/api/ribctl/__wip/data/msa_classes_proteovision/{}/{}_ribovision.fasta'.format(subunit, polymer_class)

path      = datapath('LSU',nomclass)
chief_msa:MSA = prd.parseMSA(path)

protein_id = "CAL18894.1"
url        = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id={protein_id}&retmode=json"

seq:Sequence
labels = []
seqs   = []
for seq in chief_msa:
    nih_id = seq.getLabel().split('|')[1]
    sequence_proper = barr2str(seq.getArray())
    seqs.append(sequence_proper)
    labels.append(nih_id)
    # print("label set to ", seq.getLabel())
    # print(seq)
    # response        = requests.get(url)
    # if response.status_code == 200:
    #     protein_summary = response.json()
    #     taxid           = protein_summary["result"][protein_id]["taxid"]
    # else:
    #     print("Error: ", response.status_code)


print("edited labels")
prd.writeMSA(path, MSA(np.array(seqs),nomclass,labels=labels))

# response = requests.get(url)
# protein_summary = response.json()
# taxid = protein_summary["result"][protein_id]["taxid"]
# print(taxid)

# for seq in chief_msa:
#     print(seq)
