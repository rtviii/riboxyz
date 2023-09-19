from functools import reduce
from io import StringIO
from itertools import tee
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
from ribctl.lib.classification import classify_sequence
from ribctl.lib.msalib import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.lib.ribosome_types.types_ribosome import ProteinClass, ProteinClassEnum
from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
hmm_cachedir = ASSETS['__hmm_cache']


import sys
print ("Number of arguments:", len(sys.argv), "arguments")
print ("Argument List:", str(sys.argv))


def process_struct(rcsb_id:str):
    rib            = RibosomeAssets(rcsb_id).profile()
    organism_taxid = rib.src_organism_ids[0]
    prots          = rib.proteins

    for rp in prots:
        print("processing protein {} nomenclature {} | {}".format( rp.auth_asym_id,rp.nomenclature, rp.rcsb_pdbx_description,))
        seq = rp.entity_poly_seq_one_letter_code_can
        x   = classify_sequence(seq, organism_taxid, ProteinClassEnum)
        if x not in rp.nomenclature and len(rp.nomenclature) != 0:
            print(">>>> Discovered incongruent nomenclature for protein {} nomenclature {} classification {}".format(rp.rcsb_pdbx_description, rp.nomenclature, x))
        if len(rp.nomenclature) == 0:
            print("{}".format(x))

if sys.argv[1] =="s":
    process_struct(sys.argv[2].upper())
else:
    for rcsb_id in RibosomeAssets.list_all_structs()[:10]:
        process_struct(rcsb_id)


# for candidate_class in list_ProteinClass:
#     fasta_path   = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{candidate_class}.fasta")
#     records      = Fasta(fasta_path)
#     ids          = records.all_taxids()
#     phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), n_neighbors=10)
#     seqs         = records.pick_taxids(phylo_nbhd)
#     seqs_aligned = muscle_align_N_seq( iter(seqs))
#     seqs_aligned1, seqs_aligned2 = tee(seqs_aligned)    

#     seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs_aligned1]

#     cached_name = "class_{}_taxid_{}.hmm".format(candidate_class, organism_taxid)
#     alphabet      = pyhmmer.easel.Alphabet.amino() 
#     # alphabet      = pyhmmer.easel.Alphabet.rna() 
#     builder       = pyhmmer.plan7.Builder(alphabet)
#     background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
#     anonymous_msa = pyhmmer.easel.TextMSA(bytes(cached_name, 'utf-8'),sequences=seq_tuples)
#     hmm, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)
    
#     if not os.path.isfile(os.path.join(hmm_cachedir, cached_name)):
#         with open(os.path.join(hmm_cachedir, cached_name), "wb") as hmm_file:
#             hmm.write(hmm_file)
#             print("Wrote `{}` to `{}`".format(cached_name, hmm_cachedir))
#     else:
#         print(os.path.join(hmm_cachedir, cached_name) + " exists. ")
    


#? extract sequences 
#? store (create an asset class in __init__), possibly compress
#? create method to pick phyl. nbhd, create an MSA, create an HMM
#? ->evaluate as before