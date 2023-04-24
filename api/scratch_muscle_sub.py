import os
import subprocess
from Bio import SeqIO
from Bio import SeqRecord
from prody import MSA, Sequence, MSAFile,parseMSA

from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.msa.msalib import msa_class_proteovision_path, seq_to_fasta


rcsb_id = "5AFI"
p = RibosomeAssets("5AFI").get_chain_by_polymer_class("bL25")


bl25 = p.entity_poly_seq_one_letter_code_can
seq_to_fasta('5afi', bl25, '5afi_bl25.fasta')


# ./api/ribctl/muscle3.8 -profile -in1 ./api/ribctl/assets/protein_classes_msa_proteovision/LSU/bL25_ribovision.fasta -in2 5afi_bl25.fasta -quiet
class_profile = msa_class_proteovision_path('bL25')
print(class_profile)
target_seq = '/home/rxz/dev/docker_ribxz/5afi_bl25.fasta'

from Bio.Seq import Seq
_seq          = bl25.replace("\n", "")
seq_record    = SeqRecord.SeqRecord(Seq(_seq).upper())
seq_record.id = seq_record.description = rcsb_id
bl25 = seq_record.format('fasta')

cmd = [
    '/home/rxz/dev/docker_ribxz/api/ribctl/muscle3.8',
    '-profile',
    '-in1',
    class_profile,
    '-in2',
    # target_seq,
    '-',
    '-quiet']

process = subprocess.Popen(cmd,
                           stdout=subprocess.PIPE,
                           stdin=subprocess.PIPE,
                           stderr=subprocess.PIPE, env=os.environ.copy())

stdout, stderr = process.communicate(input=bl25.encode())

print(stdout.decode())
print(stderr.decode())

# process.wait()
# print(stdout)
