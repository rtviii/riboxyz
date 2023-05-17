import os
from pyhmmer import hmmsearch,utils
import pyhmmer
from pyhmmer.easel import SequenceFile, Sequence,DigitalSequence, Alphabet

from api.ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from api.scratch_tunnel_workflow import BACTERIAL

HMM_PROFILES = "/home/rxz/dev/docker_ribxz/api/prot_hmms/"
bl12path     = "/home/rxz/dev/docker_ribxz/api/3J7Z_bL12.fasta"
# class_hmm    = "/home/rxz/dev/docker_ribxz/api/prot_hmms/bL12.hmm"


def compare_seq_to_hmm(seqpath, hmmpath):
    with pyhmmer.plan7.HMMFile(hmmpath) as hmm_file:
        hmm = hmm_file.read()

    with pyhmmer.easel.SequenceFile(seqpath, digital=True, alphabet=Alphabet.amino()) as seq_file:
        sequences = seq_file.read_block()

    pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
    return [ hit for hit in  pipeline.search_hmm(hmm, sequences) ]

for rcsb_id in BACTERIAL:
    


for i in list_ProteinClass:
    class_hmm = os.path.join(HMM_PROFILES, f"{i}.hmm")
    compare_seq_to_hmm(bl12path, class_hmm)




# seq          = SequenceFile("/home/rxz/dev/docker_ribxz/api/3J7Z_bL12.fasta", digital=True)

# sequence_file = SequenceFile(sequence_file_path)
# sequence = sequence_file.read_one_sequence()

# for hits in search:
#     print(f"HMM {hits.E} found {len(hits)} hits in the target sequences")

# alphabet      = pyhmmer.easel.Alphabet.amino()
# builder       = pyhmmer.plan7.Builder(alphabet)
# background    = pyhmmer.plan7.Background(alphabet)
# anonymous_msa = pyhmmer.easel.TextMSA(sequences=[])

# hmm, _, _ = builder.build_msa(msa, background)


# pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)
# with pyhmmer.easel.SequenceFile("data/seqs/LuxC.faa", digital=True, alphabet=alphabet) as seq_file:
#     hits = pipeline.search_hmm(hmm, seq_file)