from io import StringIO
from itertools import tee
import os
import subprocess
from tempfile import NamedTemporaryFile
from typing import Generic, Iterator, Literal
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.msalib import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import RNAClass, list_ProteinClass
from ribctl.lib.ribosome_types.types_ribosome import RNA, Protein, ProteinClass
from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
hmm_cachedir = ASSETS['__hmm_cache']



def rp_hmm_dict_init() ->dict[ProteinClass, HMM]: 
    "Retrieve dictionary of HMMs for each protein class (to compare an incoming seq against each hmm)"
    prot_hmms_dict = {}
    for hmm_class in list_ProteinClass:
        class_hmm = os.path.join(ASSETS["hmm_ribosomal_proteins"], f"{hmm_class}.hmm")
        with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
            hmm                       = hmm_file.read()
            prot_hmms_dict[hmm_class] = hmm
    return prot_hmms_dict

def seq_prot_hmm_evalue(seq:str, hmm:HMM):
    """Fit a sequence to a given HMM"""
    AMINO = Alphabet.amino()
    seq_  = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
    dsb   = DigitalSequenceBlock(AMINO, [seq_.digitize(AMINO)])
    return pyhmmer.plan7.Pipeline(alphabet=AMINO).search_hmm(hmm,dsb)

def seq_prot_against_protclasses(seq:str, hmm_dict:dict)->dict[ProteinClass, list[float]]:
    """Fit a sequence against all protein classes simultaneously"""
    _ = {}
    for (prot_class, hmm) in hmm_dict.items():
        result = seq_prot_hmm_evalue(seq, hmm)
        _.update({prot_class: [] if len(result) == 0 else list(map(lambda x: x.evalue, result))})
    return _


def classify_sequence(seq:str, organism:int, candidate_category: Literal["ribosomal_protein", "ribosomal_rna", "factor","trna"] ):
    

    ...

def classify_ribosomal_rna():
    ...

def classify_factor():
    ...

def classify_trna():
    ...

rib            = RibosomeAssets('3J7Z').profile()
organism_taxid = rib.src_organism_ids[0]
prots          = RibosomeAssets('3J7Z').profile().proteins


print(prots[14])


classify_ribosomal_protein(prots[14].entity_poly_seq_one_letter_code_can, organism_taxid)