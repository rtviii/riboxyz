# TODO

# - grab a struct->(taxid)
# - it's rrna
# - construct its phylobnhd msa
# - run through the hmm


from pprint import pprint
import pyhmmer
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libhmm import HMMs, hmm_create, seq_evaluate_v_hmm
from ete3 import NCBITaxa

from ribctl.lib.libmsa import Fasta
ALPHABETRNA = pyhmmer.easel.Alphabet.rna()

ncbi    = NCBITaxa()

RCSB_ID = "4W29"
ra      = RibosomeAssets(RCSB_ID)
taxid   = ra.get_taxids()[0][0]
chain   = ra.get_chain_by_polymer_class("23SrRNA")
seq     = chain.entity_poly_seq_one_letter_code_can

seq_record = chain.to_SeqRecord()


afasta =Fasta('23SrRNA_neighbors_of_4W29.fasta')

hmm =hmm_create('23SrRNA_neighbors_of_4W29', afasta.records, ALPHABETRNA)

pprint(seq)
query_seq = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq).digitize(ALPHABETRNA)
# hit       = seq_evaluate_v_hmm(query_seq,ALPHABETRNA, hmm)[0]


s = pyhmmer.hmmscan([ query_seq ],[ hmm ], background =pyhmmer.plan7.Background(ALPHABETRNA))
print(s)