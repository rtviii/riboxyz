import os
from pyhmmer import hmmsearch,utils
import pyhmmer
from pyhmmer.easel import SequenceFile, Sequence,DigitalSequence, Alphabet

HMM_PROFILES = "./prot_hmms"


seq       = "/home/rxz/dev/docker_ribxz/api/3J7Z_bL12.fasta"
class_hmm = "/home/rxz/dev/docker_ribxz/api/prot_hmms/bL12.hmm"


sequence_file = SequenceFile(seq)
sequence = sequence_file.read()
print(sequence)

with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
    hmm = hmm_file.read()

with pyhmmer.easel.SequenceFile(seq, digital=True, alphabet=Alphabet.amino()) as seq_file:
    sequences = seq_file.read_block()

pipeline = pyhmmer.plan7.Pipeline(hmm.alphabet)
hits     = pipeline.search_hmm(hmm, sequences)


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