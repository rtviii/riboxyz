import os
from pprint import pprint
import subprocess
from Bio import SeqRecord
from prody import Polymer
from Bio.Seq import Seq
from Bio import AlignIO
from io import StringIO
import prody
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import RNA, PolymericFactor, Protein, ProteinClass
from api.ribctl.msa.msalib import msa_class_proteovision_path, seq_to_fasta



def get_fasta_string(chain:RNA|Protein| PolymericFactor)->str:
    fasta_description      = "[{}.{}] {} |{}| {}".format(chain.parent_rcsb_id,chain.auth_asym_id, chain.src_organism_names[0], "",  chain.src_organism_ids[0])
    _seq                   = chain.entity_poly_seq_one_letter_code_can.replace("\n","")
    seq_record             = SeqRecord.SeqRecord(Seq(_seq).upper())
    seq_record.id          = fasta_description
    seq_record.description = ""
    return seq_record.format('fasta')


# ./api/ribctl/muscle3.8 -profile -in1 ./api/ribctl/assets/protein_classes_msa_proteovision/LSU/bL25_ribovision.fasta -in2 5afi_bl25.fasta -quiet
# target_seq = '/home/rxz/dev/docker_ribxz/5afi_bl25.fasta'
# print(class_profile)


rcsb_id    :str         = "5AFI"
poly_class:ProteinClass = "bL25"

def prot_class_msa_extend(rcsb_id:str, poly_class:ProteinClass)->AlignIO.MultipleSeqAlignment:
    R                 = RibosomeAssets(rcsb_id)
    chain             = R.get_chain_by_polymer_class(poly_class)
    if chain is None:
        raise LookupError("Could not find chain in {} for protein class: {}".format(rcsb_id,poly_class))

    fasta_target  = get_fasta_string(chain)
    class_profile = msa_class_proteovision_path(poly_class)

    cmd = [
        '/home/rxz/dev/docker_ribxz/api/ribctl/muscle3.8',
        '-profile',
        '-in1',
        class_profile,
        '-in2',
        '-',
        '-quiet']

    process = subprocess.Popen(cmd,
                               stdout=subprocess.PIPE,
                               stdin=subprocess.PIPE,
                               stderr=subprocess.PIPE, env=os.environ.copy())

    stdout, stderr = process.communicate(input=fasta_target.encode())
    out,err = stdout.decode(),stderr.decode()
    process.wait()

    msa_file = StringIO(out)
    msa      = AlignIO.read(msa_file, "fasta")

    return msa
msa = prot_class_msa_extend(rcsb_id,poly_class)

for s in msa:
    print(s.seq)