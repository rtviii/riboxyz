import argparse
import os
import subprocess
import sys
import prody as prd
from prody import calcShannonEntropy

from api.ribctl.msa.msalib import muscle_combine_profile

# def muscle_combine_profiles(msa_path1: str, msa_path2: str, out_filepath: str):
#     """Combine two MSA-profiles into a single one. Used here to "append" a target sequence two the ribovision alignment. """
#     cmd = ['/home/rxz/dev/docker_ribxz/cli/scripts/muscle3.8', '-profile','-in1', msa_path1, '-in2', msa_path2, '-out', out_filepath]
#     subprocess.Popen(cmd,
#                       stdout=subprocess.PIPE,
#                       stderr=subprocess.PIPE,
#                       env=os.environ.copy()
#                       ).wait()
#     sys.stdout.flush()

def barr2str (bArr):
    return ''.join([ x.decode("utf-8") for x in bArr])
    
def compare_entropies(msa_path1: str, msa_path2: str):
    msa_file     = prd.parseMSA(msa_path1)
    msa_file_ctl = prd.parseMSA(msa_path2)
    H            = sum(calcShannonEntropy(msa_file, omitgaps=False))
    H_w_ctl      = sum(calcShannonEntropy(msa_file_ctl, omitgaps=False))
    return  H, H_w_ctl

prd.confProDy(verbosity='none')

chief_msa_path = 'PV_uS5.fasta'
chief_msa      = prd.parseMSA(chief_msa_path)

suspects = [ 
'bL34',
'bS2',
'uL2',
'uS3',
'uS5'
 ]

paths = [*map(lambda _: "5afi_unclassified_{}.fasta".format(_),suspects)]

for control in suspects:
    control_path = (lambda _: "5afi_unclassified_{}.fasta".format(_))(control)
    muscle_combine_profile(chief_msa_path, control_path, '{}with_{}.fasta'.format(chief_msa_path,control))

for control in suspects:
    print("---------")
    # print(f"Comparing 'PV_uL2.fasta' against 'PV_uL2_with_{control}.fasta'") 
    print(f"Comparing 'PV_uS5.fasta' against 'PV_uS5_with_{control}.fasta'") 
    h, hwctl = compare_entropies(chief_msa_path, '{}with_{}.fasta'.format(chief_msa_path,control))
    print("Unpetrurbed entropy: {}".format(h))
    print("Perturbed entropy: {}".format(hwctl))