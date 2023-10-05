import json
import os
import pickle
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Generic, Iterator, Literal, LiteralString, Optional, Tuple, TypeVar
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
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock, DigitalSequence
from pyhmmer.plan7 import Pipeline, HMM 
import logging
import concurrent.futures
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

#? Constructon


def fasta_phylogenetic_correction(candidate_class:ProteinClass|RNAClass|LifecycleFactorClass, organism_taxid:int, max_n_neighbors:int=10)->Iterator[SeqRecord]:

    """Given a candidate class and an organism taxid, retrieve the corresponding fasta file, and perform phylogenetic correction on it."""

    if candidate_class in ProteinClass:
        fasta_path = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")

    elif candidate_class in RNAClass:
        fasta_path = os.path.join(ASSETS["fasta_rna"], f"{candidate_class.value}.fasta")

    elif candidate_class in LifecycleFactorClass:
        if candidate_class in ElongationFactorClass:
            factor_class_path = ASSETS["fasta_factors_elongation"]
        elif candidate_class in InitiationFactorClass:
            factor_class_path = ASSETS["fasta_factors_initiation"]
        else:
            raise Exception("Phylogenetic correction: Unimplemented factor class")

        fasta_path = os.path.join(factor_class_path, f"{candidate_class.value}.fasta")
    else:
        raise Exception("Phylogenetic correction: Unimplemented candidate class")

    records      = Fasta(fasta_path)
    ids          = records.all_taxids()
    phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), max_n_neighbors)
    seqs         = records.pick_taxids(phylo_nbhd)
    return iter(seqs)

def seq_evaluate_v_hmm(seq:str,alphabet:Alphabet, hmm:HMM):
    """Fit a sequence to a given HMM"""
    seq_  = pyhmmer.easel.TextSequence(name=b"template", sequence=seq)
    dsb   = DigitalSequenceBlock(alphabet, [seq_.digitize(alphabet)])
    return pyhmmer.plan7.Pipeline(alphabet=alphabet, T=100).search_hmm(hmm,dsb)

def seq_evaluate_v_hmm_dict(seq:str,alphabet:Alphabet, hmm_dict:dict)->dict[PolymerClass, list[float]]:
    """Fit a sequence against all protein classes simultaneously"""
    _ = {}
    for (candidate_class, hmm) in hmm_dict.items():
        result = seq_evaluate_v_hmm(seq,alphabet, hmm)
        _.update({candidate_class: [] if len(result) == 0 else list(map(lambda x: { "evalue":x.evalue, "score":x.score }, result))})
    return _


def pick_best_hmm_hit(matches_dict:dict[PolymerClass, list[float]], chain_info:Polymer)->PolymerClass | None:
    """Given a dictionary of sequence-HMMe e-values, pick the best candidate class"""
    results = []
    if len([ item for x in list(matches_dict.values()) for item in x ]) > 0:
        pprint(matches_dict)
    for (candidate_class, match) in matches_dict.items():
        if len(match) == 0:
            continue
        else:
            results.append({"candidate_class":candidate_class,"match": match})

    if len(results) == 0 :
        # logging.warning("Chain {}.{} : Did not match any of the candidate HMM models.".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id))
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

    elif candidate_category == RNAClass:
        candidates_dict = candidates_dict if candidates_dict is not None else hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results         = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.rna(), candidates_dict)
        return results

    elif candidate_category == LifecycleFactorClass:
        candidates_dict = candidates_dict if candidates_dict is not None else hmm_dict_init__candidates_per_organism(candidate_category, organism)
        results         = seq_evaluate_v_hmm_dict(seq, pyhmmer.easel.Alphabet.amino(), candidates_dict)
        return results
    else:
        raise Exception("Classify sequence: Unimplemented candidate category")

#! Implementations ------------------------------

def classify_subchain(chain: typing.Union[Protein, RNA, LifecycleFactor] , candidates_dict:dict[PolymerClass, HMM]|None=None)->Tuple[str, PolymerClass|None]:
    # logging.debug("Task for chain {}.{} (Old nomenclature {}) | taxid {}".format( chain.parent_rcsb_id,chain.auth_asym_id, chain.nomenclature, chain.src_organism_ids[0]))
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

def hmm_cache(hmm:HMM):
    name     = hmm.name.decode('utf-8')
    filename = os.path.join(hmm_cachedir, name)

    if not os.path.isfile(filename):
        with open(filename, "wb") as hmm_file:
            hmm.write(hmm_file)
            print("Wrote `{}` to `{}`".format(filename, hmm_cachedir))
    else:
        ...

def hmm_check_cache(candidate_class: ProteinClass | RNAClass | LifecycleFactorClass, organism_taxid:int)->HMM | None:
    hmm_path = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)
    if os.path.isfile(os.path.join(hmm_cachedir, hmm_path)):
        hmm_path = os.path.join(hmm_cachedir, hmm_path)
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            HMM = hmm_file.read()
            return HMM
    else:
        return None

def hmm_dict_init__candidates_per_organism(candidate_classes: list[ PolymerClass ] ,organism_taxid:int)->dict[PolymerClass, HMM]:

    _ ={}
    for candidate_class in candidate_classes:
        if candidate_class in ProteinClass:
            _.update({ candidate_class.value: HMMScanner.hmm_produce(candidate_class, organism_taxid) })

        elif candidate_class in RNAClass:
            _.update({ candidate_class.value: HMMScanner.hmm_produce(candidate_class, organism_taxid) })

        elif candidate_class in LifecycleFactorClass:
            _.update({ candidate_class.value: HMMScanner.hmm_produce(candidate_class, organism_taxid) })
        else:
            raise Exception("hmm_dict_init__candidates_per_organism: Unexpected candidate class: {}".format(candidate_class))
    return {}

def hmm_create(name:str, seqs:Iterator[SeqRecord], alphabet:Alphabet)->HMM:
    """Create an HMM from a list of sequences"""

    seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs]
    builder       = pyhmmer.plan7.Builder(alphabet)
    background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
    anonymous_msa = pyhmmer.easel.TextMSA(bytes(name, 'utf-8'),sequences=seq_tuples)
    HMM, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)

    return HMM

def hmm_produce(candidate_class: ProteinClass | RNAClass | LifecycleFactorClass, organism_taxid:int, no_cache:bool=False, max_seed_seq:int=10)->HMM:  # type: ignore
    """Produce an organism-specific HMM. Retrieve from cache if exists, otherwise generate and cache."""
    if ( hmm := hmm_check_cache(candidate_class, organism_taxid) ) != None and not no_cache:
        return hmm
    else:
        if candidate_class in ProteinClass or candidate_class in LifecycleFactorClass:
            alphabet = pyhmmer.easel.Alphabet.amino()
        elif candidate_class in LifecycleFactorClass: 
            alphabet = pyhmmer.easel.Alphabet.rna()
        else:
            raise Exception("hmm_produce: Unimplemented candidate class")
                
        seqs   = fasta_phylogenetic_correction(candidate_class, organism_taxid, max_n_neighbors=max_seed_seq)
        seqs_a = muscle_align_N_seq(iter(seqs))
        name   = "{}".format(candidate_class.value)
        HMM    = hmm_create(name, seqs_a, alphabet)

        if not no_cache:
            hmm_cache(HMM)
        return HMM

class HMMScanner():
    # https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.HMM
    """
    STEPS (bottom-up):
    - given a taxonomic id, produce a set of HMMs 
    - each HMM is representative of a polymer class parametrized by a seed MSA 
    - a sequence is searched against all HMMs in the registry
    """

    seed_sequences : dict[PolymerClass, list[SeqRecord]] = {}
    hmms_registry  : dict[PolymerClass, HMM]             = {}


    def __init__(self, tax_id:int, candidate_classes:list[PolymerClass],  no_cache:bool=False, max_seed_seqs:int=10) -> None:

        self.organism_tax_id = tax_id

        def __load_seed_sequences(candidate_class:PolymerClass,tax_id:int, max_seed_seqs:int):
            """This is just for housekeeping. Keeping sequences with which the HMMs were seeded as state on the classifier."""
            return [*fasta_phylogenetic_correction(candidate_class, tax_id, max_n_neighbors=max_seed_seqs)]

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

        for candidate in candidate_classes:



# !-
        for candidate_class in candidate_classes:
            seqs                                       = [*fasta_phylogenetic_correction(candidate_class, tax_id, max_n_neighbors=max_seed_seqs)]
            self.seed_sequences[candidate_class.value] = seqs
            self.hmms_registry[candidate_class.value]  = hmm_produce(candidate_class, tax_id, no_cache=no_cache, max_seed_seq=max_seed_seqs)
            print("Loded HMM: {}".format(candidate_class))
# !-

    def scan(self, alphabet:pyhmmer.easel.Alphabet, target_seqs:list[SeqRecord]):
        """Construct a scan pipeline for the current classifier"""

        # convert seq records to a list of pyhmmer.DigitalSequences
        query_seqs = []
        for seq_record in target_seqs:
            query_seq  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
            query_seqs.append(query_seq.digitize(alphabet))
        scans = list(pyhmmer.hmmscan(query_seqs,[*self.hmms_registry.values()] ))
        return scans

    def info(self)->dict:
        """Get info for the constituent HMMs in the current classifier"""
        _ = {
             "seed_seqs": {
                str(k.value): [{"seq":str(seqrecord.seq), "tax_id":seqrecord.id} for seqrecord in v] for (k, v) in self.seed_sequences.items()
             },
             "hmms_registry": {
                    str( k.value ): {"name": v.name.decode(),  "M":v.M, "nseq":v.nseq, "nseq_effective":v.nseq_effective}
                    for (k, v) in self.hmms_registry.items()
                },
             }
        return _


class HMMClassifier():

    organism_scanners : dict[int, HMMScanner] = {}
    chains            : list[Polymer] = []
    alphabet          : pyhmmer.easel.Alphabet
    candidate_classes : list[PolymerClass]

    bitscore_threshold: int = 35 # TODO: unused as of now
    report      : dict

    def __init__(self, chains: list[Polymer], alphabet:pyhmmer.easel.Alphabet) -> None:
        self.chains         = chains
        self.alphabet       = alphabet
        self.report         = {}

        if alphabet == Alphabet.amino():
            self.candidate_classes  = [pc for pc in [*list(ProteinClass), *list(LifecycleFactorClass)]]

        elif alphabet == Alphabet.rna():
            self.candidate_classes  = [pc for pc in [*list(RNAClass)]]

    def pick_best_hit(self, hits:list[dict]) -> list[PolymerClass]:
        # TODO: look at filtering allunder self.bitscore threshold (while scaling the scores by the [average of the lowest half * 1.5]/[average of the lowest half])
        if len(hits) == 0:
            return []

        sorted_by_biscore = sorted(hits, key=lambda x: x["bitscore"], reverse=True)
        return sorted_by_biscore[0]['class_name'] if sorted_by_biscore[0]['bitscore'] > self.bitscore_threshold else []


        # return [ sorted(hits, key=lambda x: x["hit.score"])[0]["hit.name"].decode() ] if len(hits) >0 else []

    def scan_chains(self)->None:

        for chain in self.chains:
            organism_taxid = chain.src_organism_ids[0]

            # -- pick scanner'|
            if organism_taxid not in self.organism_scanners:
                    hmmscanner = HMMScanner( organism_taxid, self.candidate_classes, no_cache     = True, max_seed_seq = 5)
                    print("Added scanner {} for organism {}".format( hmmscanner.__str__(), organism_taxid ))
                    self.organism_scanners[organism_taxid] = hmmscanner

            else:
                hmmscanner = self.organism_scanners[organism_taxid]
            # -- pick scanner.|
 

            self.report[chain.auth_asym_id] = []
            # -- convert seq to easel format
            seq_record = chain.to_SeqRecord()
            query_seq  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
            query_seqs = [query_seq.digitize(self.alphabet)]
            # -- convert seq to easel format
            
            for scan in list(pyhmmer.hmmscan(query_seqs,[*hmmscanner.hmms_registry.values()])):
                for hit in [*scan]:
                   d_hit = {
                    "fasta_seed"        : [str(seqrecord.seq) for seqrecord in hmmscanner.seed_sequences[hit.name.decode()]],
                    "seed_organism_ids" : [seqrecord.id for seqrecord in hmmscanner.seed_sequences[hit.name.decode()]],
                    "class_name"        : hit.name.decode(),                                                             # <- comes from the emitting hmm
                    "evalue"            : hit.evalue,
                    "bitscore"          : hit.score,
                    "target_organism_id": organism_taxid,
                    "consensus"         : hmmscanner.hmms_registry[hit.name.decode()].consensus,
                    "domains"           : [( d.score, d.c_evalue, d.env_from, d.env_to ) for d in hit.domains]
                    }

                   self.report[chain.auth_asym_id].append(d_hit)

    def produce_classification(self)->dict[str,list]:
        classes = {}
        for (auth_asym_id,hits) in self.report.items():
            classes[auth_asym_id] = self.pick_best_hit(hits)
        return classes
   
    def write_classification_report(self,write_path:str)->None:
        with open(write_path, "w") as outfile:
            json.dump(self.report, outfile)


