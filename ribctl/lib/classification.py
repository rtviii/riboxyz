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
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import PolymericFactorClass, RNAClass, list_ProteinClass
from ribctl.lib.ribosome_types.types_ribosome import RNA, Polymer, PolymerClass, PolymerClass_, PolymericFactor, Protein, ProteinClass, ProteinClassEnum, RNAClassEnum
# from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM 
import logging
import concurrent.futures

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
        for rc in RNAClassEnum:
            _.update({ rc.value: hmm_produce(rc, organism_taxid) })

    elif candidate_category == PolymericFactorClass:
        raise Exception("Not implemented yet: PolymericFactorClass hmmdictinit")
    else:
        return {}
    
    return _

def pick_best_class(matches_dict:dict[PolymerClass_, list[float]], chain_info:Polymer)->PolymerClass_ | None:
    """Given a dictionary of sequence-HMMe e-values, pick the best candidate class"""
    results = []
    for (candidate_class, matches) in matches_dict.items():
        if len(matches) == 0:
            continue
        else:
            results.append((candidate_class, matches))
    if len(results) == 0 :
        logging.warning("Chain {}.{} : Did not match any of the candidate HMM models.".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id))
        return None
    if len(results) > 1 :
        # if more than 1 match, pick the smallest and ring alarms if the next smallest is within an order of magnitude.
        results = sorted(results, key=lambda match_kv: match_kv[1])
        if abs(math.log10(results[0][1][0]/results[1][1][0])) < 2:
            logging.warning("{}.{} : Multiple sensible matches detected, picked {} (smallest e-value) : \n {}".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id,results[0][0], results))
    return results[0][0]

def classify_sequence(seq:str, organism:int, candidate_category:typing.Union[RNAClassEnum, ProteinClassEnum], candidates_dict:dict[PolymerClass_, HMM]|None=None)->dict[PolymerClass_, list[float]]:

    if candidate_category == ProteinClassEnum:
        candidates_dict = candidates_dict if candidates_dict is not None else hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results         = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.amino(), candidates_dict)
        return results

    if candidate_category == RNAClassEnum:
        #TODO
        raise Exception("Not implemented yet")

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

            seqs = fasta_phylogenetic_correction(candidate_class, organism_taxid, n_neighbors=10)
            seqs_a = muscle_align_N_seq(iter(seqs))

            cached_name = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)
            alphabet    = pyhmmer.easel.Alphabet.amino()
            HMM         = hmm_create(cached_name, seqs_a, alphabet)

            hmm_cache(HMM)
            
            return HMM

        if candidate_class in RNAClassEnum:
            raise Exception("Not implemented yet")
            ...

#! Implementations ------------------------------

def classify_subchain(chain: typing.Union[Protein, RNA, PolymericFactor] , candidates_dict:dict[PolymerClass_, HMM]|None=None)->Tuple[str, PolymerClass_|None]:
    logging.debug("Task for chain {}.{} (Old nomenclature {}) | taxid {}".format( chain.parent_rcsb_id,chain.auth_asym_id, chain.nomenclature, chain.src_organism_ids[0]))
    if type(chain) == RNA:
        assigned = classify_sequence(chain.entity_poly_seq_one_letter_code_can, chain.src_organism_ids[0], RNAClassEnum, candidates_dict=candidates_dict)

    elif type(chain) == Protein:
        assigned = classify_sequence(chain.entity_poly_seq_one_letter_code_can, chain.src_organism_ids[0], ProteinClassEnum, candidates_dict=candidates_dict)

    elif type(chain)== PolymericFactor:
        return (chain.auth_asym_id, None)
    else:
        raise Exception("Invalid chain type")
    
    assigned = pick_best_class(assigned, chain)
    return (chain.auth_asym_id, assigned)

# def classify_subchains(targets:list[Protein]|list[RNA]|list[PolymericFactor])->dict[str, PolymerClass_]:
def classify_subchains(targets:list[typing.Union[Protein,RNA,PolymericFactor]],candidate_category:PolymerClass_)->dict[str, PolymerClass_]:
    # This is a dict of dicts (a "registry", sigh) to only do i/o once per organism.(Pass it down)
    # It's a dictionary because some chains within a structure might originate in different organisms. 
    organisms_registry: dict[int,dict[PolymerClass_, HMM]]= {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
        tasks   = []
        results = []
        for chain in targets:
            if chain.src_organism_ids[0] not in [*organisms_registry.keys()]:
                organisms_registry[chain.src_organism_ids[0]] = hmm_dict_init__candidates_per_organism(candidate_category, chain.src_organism_ids[0])

            future = executor.submit(classify_subchain, chain, organisms_registry[chain.src_organism_ids[0]])
            tasks.append(future)

        for future in concurrent.futures.as_completed(tasks):
            results.append(future.result())

        return {k:v for (k,v) in results}



