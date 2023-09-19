from enum import Enum
from io import StringIO
from itertools import tee
import math
import os
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Generic, Iterator, Literal, LiteralString, TypeVar
import typing
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.msalib import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import RNAClass, list_ProteinClass
from ribctl.lib.ribosome_types.types_ribosome import RNA, PolymerClass, PolymerClass_, Protein, ProteinClass, ProteinClassEnum, RNAClassEnum
# from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 

hmm_cachedir = ASSETS['__hmm_cache']



# def rp_hmm_dict_init() ->dict[ProteinClass, HMM]: 
    # "Retrieve dictionary of HMMs for each protein class (to compare an incoming seq against each hmm)"
    # prot_hmms_dict = {}
    # for hmm_class in list_ProteinClass:
    #     class_hmm = os.path.join(ASSETS["hmm_ribosomal_proteins"], f"{hmm_class}.hmm")
    #     with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
    #         hmm                       = hmm_file.read()
    #         prot_hmms_dict[hmm_class] = hmm
    # return prot_hmms_dict




    # fasta_path   = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{candidate_class}.fasta")
    # records      = Fasta(fasta_path)
    # ids          = records.all_taxids()
    # phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), n_neighbors=10)
    # seqs         = records.pick_taxids(phylo_nbhd)
    # seqs_aligned = muscle_align_N_seq( iter(seqs))
    # seqs_aligned1, seqs_aligned2 = tee(seqs_aligned)    

    # seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs_aligned1]



# def load_hmm_from_cache(path:str)->HMM:
#     if os.path.isfile(path):
#         with pyhmmer.plan7.HMMFile(path) as hmm_file:
#             hmm = hmm_file.read()
#     else:
        # generate the given hmm


def rp_hmm_dict_init(organim_taxid:int) ->dict[ProteinClass, HMM]: 
    "Retrieve dictionary of HMMs for each protein class (to compare an incoming seq against each hmm)"
    prot_hmms_dict = {}
    for hmm_class in list_ProteinClass:
        class_hmm = os.path.join(ASSETS["hmm_ribosomal_proteins"], f"{hmm_class}.hmm")
        with pyhmmer.plan7.HMMFile(class_hmm) as hmm_file:
            hmm                       = hmm_file.read()
            prot_hmms_dict[hmm_class] = hmm
    return prot_hmms_dict

#!------------------ Former
def seq_evaluate_v_hmm(seq:str,alphabet:Alphabet, hmm:HMM):
    """Fit a sequence to a given HMM"""
    seq_  = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
    dsb   = DigitalSequenceBlock(alphabet, [seq_.digitize(alphabet)])
    return pyhmmer.plan7.Pipeline(alphabet=alphabet).search_hmm(hmm,dsb)

def seq_evaluate_v_hmm_dict(seq:str,alphabet:Alphabet, hmm_dict:dict)->dict[PolymerClass_, list[float]]:
    """Fit a sequence against all protein classes simultaneously"""
    _ = {}
    for (candidate_class, hmm) in hmm_dict.items():
        result = seq_evaluate_v_hmm(seq,alphabet, hmm)
        _.update({candidate_class: [] if len(result) == 0 else list(map(lambda x: x.evalue, result))})
    return _

def hmm_dict_init__candidates_per_organism(candidate_category:PolymerClass_,organism_taxid:int)->dict[PolymerClass_, HMM]:
    _ ={}
    if candidate_category == ProteinClassEnum:
        for pc in ProteinClassEnum:
            _.update({ pc.value: hmm_produce(pc, organism_taxid) })

    elif candidate_category == RNAClassEnum:
        raise Exception("Not implemented yet")
        for rc in RNAClassEnum:
            _.update({ rc.value: hmm_produce(rc, organism_taxid) })
    
    return _


def pick_best_candidate(matches_dict:dict[PolymerClass_, list[float]])->PolymerClass_:
    """Given a dictionary of matches, pick the best candidate class"""
    results = []
    for (candidate_class, matches) in matches_dict.items():
        if len(matches) == 0:
            continue
        else:
            results.append((candidate_class, matches))
    if len(results) == 0 :
        raise Exception("Did not detect any matches. Something went wrong.")
    if len(results) > 1 :
        # if more than 1 match, pick the smallest and ring alarms if the next smallest is within an order of magnitude.
        results = sorted(results, key=lambda match_kv: match_kv[1])
        if abs(math.log10(results[0][1]/results[1][0])) < 2:
            raise Exception("Multiple sensible matches detected. Something went wrong. \n {}".format(results))
          

    return results[0][0]



def classify_sequence(seq:str, organism:int, candidate_category:typing.Union[RNAClassEnum, ProteinClassEnum])->PolymerClass_:

    if candidate_category == ProteinClassEnum:

        candidates_dict = hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.amino(), candidates_dict)
        return pick_best_candidate(results)
    if candidate_category == RNAClassEnum:
        #TODO

        raise Exception("Not implemented yet")
        candidates_dict = hmm_dict_init__candidates_per_organism(candidate_category, organism)
        seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.rna(), candidates_dict)

def hmm_create(name:str, seqs:Iterator[SeqRecord], alphabet:Alphabet)->HMM:
    """Create an HMM from a list of sequences"""
    seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs]
    builder       = pyhmmer.plan7.Builder(alphabet)
    background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
    anonymous_msa = pyhmmer.easel.TextMSA(bytes(name, 'utf-8'),sequences=seq_tuples)
    HMM, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)
    return HMM

def fasta_phylogenetic_correction(candidate_class:ProteinClassEnum|RNAClassEnum, organism_taxid:int, n_neighbors=10)->Iterator[SeqRecord]:
    """Given a candidate class and an organism taxid, retrieve the corresponding fasta file, and perform phylogenetic correction on it."""
    if candidate_class in ProteinClassEnum:
        fasta_path = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{candidate_class.value}.fasta")
    elif candidate_class in RNAClassEnum:
        fasta_path = os.path.join(ASSETS["fasta_ribosomal_rna"], f"{candidate_class.value}.fasta")
    else:
        raise Exception("Invalid candidate class")

    records      = Fasta(fasta_path)
    ids          = records.all_taxids()
    phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), n_neighbors)
    seqs         = records.pick_taxids(phylo_nbhd)
    return iter(seqs)

def hmm_cache(hmm:HMM):
    name     = hmm.name.decode('utf-8')
    filename = os.path.join(hmm_cachedir, name)

    if not os.path.isfile(filename):
        with open(filename, "wb") as hmm_file:
            hmm.write(hmm_file)
            print("Wrote `{}` to `{}`".format(filename, hmm_cachedir))
    else:
        ...
        # print(filename + " exists.")

def hmm_produce(candidate_class: ProteinClassEnum | RNAClassEnum, organism_taxid:int)->HMM:  # type: ignore
    """Produce an organism-specific HMM. Retrieve from cache if exists, otherwise generate and cache."""
    hmm_path = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)

    if os.path.isfile(os.path.join(hmm_cachedir, hmm_path)):
        hmm_path = os.path.join(hmm_cachedir, hmm_path)
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            HMM = hmm_file.read()
            return HMM
    else:
        if candidate_class in ProteinClassEnum:

            # fasta_path   = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{candidate_class}.fasta")
            # records      = Fasta(fasta_path)
            # ids          = records.all_taxids()
            # phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), n_neighbors=10)
            # seqs         = records.pick_taxids(phylo_nbhd)
            # seqs_aligned = muscle_align_N_seq( iter(seqs))

            
            # seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs_aligned]
            # builder       = pyhmmer.plan7.Builder(alphabet)
            # background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
            # anonymous_msa = pyhmmer.easel.TextMSA(bytes(cached_name, 'utf-8'),sequences=seq_tuples)
            # HMM, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)

            seqs = fasta_phylogenetic_correction(candidate_class, organism_taxid, n_neighbors=10)
            seqs_a = muscle_align_N_seq(iter(seqs))

            cached_name = "class_{}_taxid_{}.hmm".format(candidate_class, organism_taxid)
            alphabet    = pyhmmer.easel.Alphabet.amino()
            HMM         = hmm_create(cached_name, seqs_a, alphabet)

            hmm_cache(HMM)
            
            #TODO: Cache it
            return HMM

        if candidate_class in RNAClassEnum:
            raise Exception("Not implemented yet")
            ...

