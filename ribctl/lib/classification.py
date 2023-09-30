from enum import Enum
from io import StringIO
from itertools import tee
import math
import os
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Generic, Iterator, Literal, LiteralString, Tuple, TypeVar
import typing
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.msalib import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_ribosome import RNA, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, Polymer, PolymerClass, LifecycleFactor, Protein, ProteinClass, ProteinClass, RNAClass
# from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
import logging
import concurrent.futures
import inspect

from ribctl.lib.util_taxonomy import taxid_domain


# Configure the logging settings
logging.basicConfig(
    level=logging.DEBUG,  # Set the logging level to DEBUG (you can adjust this)
    format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',  # Define the log message format
    handlers=[
        logging.StreamHandler(),  # Log to the console
        logging.FileHandler('classification.log')  # Log to a file named 'my_log_file.log'
    ]
)


hmm_cachedir = ASSETS['__hmm_cache']

#! Lib ------------------------------
def seq_evaluate_v_hmm(seq:str,alphabet:Alphabet, hmm:HMM):
    """Fit a sequence to a given HMM"""
    seq_  = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
    dsb   = DigitalSequenceBlock(alphabet, [seq_.digitize(alphabet)])
    return pyhmmer.plan7.Pipeline(alphabet=alphabet).search_hmm(hmm,dsb)

def seq_evaluate_v_hmm_dict(seq:str,alphabet:Alphabet, hmm_dict:dict)->dict[PolymerClass, list[float]]:
    """Fit a sequence against all protein classes simultaneously"""
    _ = {}
    for (candidate_class, hmm) in hmm_dict.items():
        result = seq_evaluate_v_hmm(seq,alphabet, hmm)
        _.update({candidate_class: [] if len(result) == 0 else list(map(lambda x: { "evalue":x.evalue, "score":x.score }, result))})
    return _

def hmm_dict_init__candidates_per_organism(candidate_category:PolymerClass,organism_taxid:int)->dict[PolymerClass, HMM]:
    _ ={}
    if candidate_category == ProteinClass:
        for pc in ProteinClass:
            _.update({ pc.value: hmm_produce(pc, organism_taxid) })

    elif candidate_category == RNAClass:
        for rc in RNAClass:
            _.update({ rc.value: hmm_produce(rc, organism_taxid) })

    elif candidate_category == LifecycleFactorClass:
        for rc in LifecycleFactorClass:
            _.update({ rc.value: hmm_produce(rc, organism_taxid) })
    else:
        return {}
    
    return _

def pick_best_hmm_hit(matches_dict:dict[PolymerClass, list[float]], chain_info:Polymer)->PolymerClass | None:
    """Given a dictionary of sequence-HMMe e-values, pick the best candidate class"""
    results = []
    for (candidate_class, match) in matches_dict.items():
        if len(match) == 0:
            continue
        else:
            results.append({"candidate_class":candidate_class,"match": match})

    if len(results) == 0 :
        logging.warning("Chain {}.{} : Did not match any of the candidate HMM models.".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id))
        return None

    if len(results) > 1 :
        # if more than 1 match, pick the smallest and ring alarms if the next smallest is within an order of magnitude.
        results = sorted(results, key=lambda match_kv: match_kv['match'][0]['score'], reverse=True)

        best_match        = results[0]['match'][0]
        second_best_match = results[1]['match'][0]
        
        if best_match == 0 and second_best_match == 0:
            logging.warning("{}.{} : Multiple sensible matches detected, picked {} (smallest e-value) : \n {}".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id,results[0][0], results))

    return results[0]['candidate_class']

def classify_sequence(seq:str, organism:int, candidate_category:typing.Union[RNAClass, ProteinClass], candidates_dict:dict[PolymerClass, HMM]|None=None)->dict[PolymerClass, list[float]]:

    if candidate_category == ProteinClass:
        candidates_dict = candidates_dict if candidates_dict is not None else hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results         = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.amino(), candidates_dict)
        return results

    if candidate_category == RNAClass:
        candidates_dict = candidates_dict if candidates_dict is not None else hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results         = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.rna(), candidates_dict)
        return results

def fasta_phylogenetic_correction(candidate_class:ProteinClass|RNAClass|LifecycleFactorClass, organism_taxid:int, max_n_neighbors=10)->Iterator[SeqRecord]:
    """Given a candidate class and an organism taxid, retrieve the corresponding fasta file, and perform phylogenetic correction on it."""

    if candidate_class in ProteinClass:
        fasta_path = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")

    elif candidate_class in RNAClass:
        fasta_path = os.path.join(ASSETS["fasta_ribosomal_rna"], f"{candidate_class.value}.fasta")

    elif candidate_class in LifecycleFactorClass:
        if candidate_class in ElongationFactorClass:
            factor_class_path = ASSETS["fasta_factors_elongation"]
        elif candidate_class in InitiationFactorClass:
            factor_class_path = ASSETS["fasta_factors_initiation"]
        else:
            raise Exception("Phylogenetic correction: Unimplemented factor class")

        domain     = taxid_domain(organism_taxid)
        fasta_path = os.path.join(factor_class_path, f"{candidate_class.value}.fasta")
    else:
        raise Exception("Phylogenetic correction: Unimplemented candidate class")


    records      = Fasta(fasta_path)
    ids          = records.all_taxids()
    phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), max_n_neighbors)
    seqs         = records.pick_taxids(phylo_nbhd)
    return iter(seqs)

def hmm_create(name:str, seqs:Iterator[SeqRecord], alphabet:Alphabet)->HMM:
    """Create an HMM from a list of sequences"""

    seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs]

    builder       = pyhmmer.plan7.Builder(alphabet)
    background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
    anonymous_msa = pyhmmer.easel.TextMSA(bytes(name, 'utf-8'),sequences=seq_tuples)

    HMM, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)
    return HMM

def hmm_cache(hmm:HMM):
    name     = hmm.name.decode('utf-8')
    filename = os.path.join(hmm_cachedir, name)

    if not os.path.isfile(filename):
        with open(filename, "wb") as hmm_file:
            hmm.write(hmm_file)
            print("Wrote `{}` to `{}`".format(filename, hmm_cachedir))
    else:
        ...

def hmm_produce(candidate_class: ProteinClass | RNAClass | LifecycleFactorClass, organism_taxid:int, no_cache:bool=False)->HMM:  # type: ignore
    """Produce an organism-specific HMM. Retrieve from cache if exists, otherwise generate and cache."""
    hmm_path = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)

    if os.path.isfile(os.path.join(hmm_cachedir, hmm_path)):
        hmm_path = os.path.join(hmm_cachedir, hmm_path)
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            HMM = hmm_file.read()
            return HMM
    else:
        if candidate_class in ProteinClass or candidate_class in LifecycleFactorClass:
            alphabet = pyhmmer.easel.Alphabet.amino()
        elif candidate_class in LifecycleFactorClass: 
            alphabet = pyhmmer.easel.Alphabet.rna()
        else:
            raise Exception("hmm_produce: Unimplemented candidate class")
                
        seqs        = fasta_phylogenetic_correction(candidate_class, organism_taxid, max_n_neighbors=10)
        seqs_a      = muscle_align_N_seq(iter(seqs))
        cached_name = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)

        HMM         = hmm_create(cached_name, seqs_a, alphabet)
        if not no_cache:
            hmm_cache(HMM)
        return HMM

#! Implementations ------------------------------

def classify_subchain(chain: typing.Union[Protein, RNA, LifecycleFactor] , candidates_dict:dict[PolymerClass, HMM]|None=None)->Tuple[str, PolymerClass|None]:
    logging.debug("Task for chain {}.{} (Old nomenclature {}) | taxid {}".format( chain.parent_rcsb_id,chain.auth_asym_id, chain.nomenclature, chain.src_organism_ids[0]))
    print("Classifying chain ", chain.auth_asym_id, chain.src_organism_ids[0])
    if type(chain) == RNA:
        assigned = classify_sequence(chain.entity_poly_seq_one_letter_code_can, chain.src_organism_ids[0], RNAClass, candidates_dict=candidates_dict)

    elif type(chain) == Protein:
        assigned = classify_sequence(chain.entity_poly_seq_one_letter_code_can, chain.src_organism_ids[0], ProteinClass, candidates_dict=candidates_dict)

    elif type(chain)== LifecycleFactor:
        assigned = classify_sequence(chain.entity_poly_seq_one_letter_code_can, chain.src_organism_ids[0], LifecycleFactorClass, candidates_dict=candidates_dict)
    else:
        raise Exception("Invalid chain type")
    
    assigned = pick_best_hmm_hit(assigned, chain)
    return (chain.auth_asym_id, assigned)

# def classify_subchains(targets:list[Protein]|list[RNA]|list[PolymericFactor])->dict[str, PolymerClass_]:
def classify_subchains(targets:list[typing.Union[Protein,RNA,LifecycleFactor]],candidate_category:PolymerClass)->dict[str, PolymerClass]:
    # This is a dict of dicts to only do IO once(on average) per organism assuming 90%+ of chains in a structure are of the same taxid.(Pass it down)
    # It's a dictionary because some chains within a structure might originate in different organisms (host system).
    hmm_organisms_registry: dict[int,dict[PolymerClass, HMM]]= {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        tasks   = []
        results = []
        for chain in targets:
            chain_organism_taxid  = chain.src_organism_ids[0]
            if chain_organism_taxid not in [*hmm_organisms_registry.keys()]:
                hmm_organisms_registry[chain_organism_taxid] = hmm_dict_init__candidates_per_organism(candidate_category, chain_organism_taxid)
            future = executor.submit(classify_subchain, chain, hmm_organisms_registry[chain_organism_taxid])
            tasks.append(future)

        for future in concurrent.futures.as_completed(tasks):
            results.append(future.result())

        return {k:v for (k,v) in results}



