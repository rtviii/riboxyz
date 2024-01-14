import json
import os
import pickle
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Generic, Iterator, Literal, Optional, Tuple, Type, TypeVar
import typing
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq 
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.libmsa import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.ribosome_types.types_ribosome import RNA, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, PolynucleotideClass, PolypeptideClass, Protein, CytosolicProteinClass, CytosolicProteinClass, tRNA
# from ribctl.etl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock, DigitalSequence
from pyhmmer.plan7 import Pipeline, HMM , TopHits
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

#? Constructon
def fasta_phylogenetic_correction(candidate_class:PolymerClass, organism_taxid:int, max_n_neighbors:int=10)->Iterator[SeqRecord]:
    """Given a candidate class and an organism taxid, retrieve the corresponding fasta file, and perform phylogenetic correction on it."""

    if candidate_class in CytosolicProteinClass:
        fasta_path = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")

    elif candidate_class in MitochondrialProteinClass:
        fasta_path = os.path.join(ASSETS["fasta_proteins_mitochondrial"], f"{candidate_class.value}.fasta")

    elif candidate_class in PolynucleotideClass:
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

    if not os.path.isfile(fasta_path):
        raise Exception("Not found: {}".format(fasta_path))

    records    = Fasta(fasta_path)
    if len( list(records.records) ) == 0 :
        raise Exception("Empty fasta file: {}".format(fasta_path))
    ids        = records.all_taxids()
    phylo_nbhd = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), max_n_neighbors)
    seqs       = records.pick_taxids(phylo_nbhd)

    return iter(seqs)

def digitize_seq_record(seq_record:SeqRecord, alphabet:Alphabet)->DigitalSequence:
    """Convert a SeqRecord to a DigitalSequence"""
    seq_  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
    return seq_.digitize(alphabet)

def seq_evaluate_v_hmm(seq:DigitalSequence,alphabet:Alphabet, hmm:HMM, T:int=35)->TopHits:
    """Fit a sequence to a given HMM"""
    dsb   = DigitalSequenceBlock(alphabet, [seq])
    return pyhmmer.plan7.Pipeline(alphabet=alphabet, T=T).search_hmm(hmm,dsb)

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

def hmm_cache(hmm:HMM):
    name     = hmm.name.decode('utf-8')
    filename = os.path.join(hmm_cachedir, name)

    if not os.path.isfile(filename):
        with open(filename, "wb") as hmm_file:
            hmm.write(hmm_file)
            print("Wrote `{}` to `{}`".format(filename, hmm_cachedir))
    else:
        ...

def hmm_check_cache(candidate_class: PolymerClass, organism_taxid:int)->HMM | None:
    hmm_path = "class_{}_taxid_{}.hmm".format(candidate_class.value, organism_taxid)
    if os.path.isfile(os.path.join(hmm_cachedir, hmm_path)):
        hmm_path = os.path.join(hmm_cachedir, hmm_path)
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            HMM = hmm_file.read()
            return HMM
    else:
        return None

def hmm_create(name:str, seqs:Iterator[SeqRecord], alphabet:Alphabet)->HMM:
    """Create an HMM from a list of sequences"""

    seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs]
    builder       = pyhmmer.plan7.Builder(alphabet)
    background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
    anonymous_msa = pyhmmer.easel.TextMSA(bytes(name, 'utf-8'),sequences=seq_tuples)
    HMM, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)

    return HMM

def hmm_produce(candidate_class: PolymerClass, organism_taxid:int, seed_sequences:list[SeqRecord], no_cache:bool=False)->Tuple[PolymerClass,HMM]:  # type: ignore
    """Produce an organism-specific HMM. Retrieve from cache if exists, otherwise generate and cache."""
    if ( hmm := hmm_check_cache(candidate_class, organism_taxid) ) != None and not no_cache:
        return (candidate_class, hmm )
    else:
        if candidate_class in CytosolicProteinClass or candidate_class in LifecycleFactorClass or candidate_class in MitochondrialProteinClass:
            alphabet = pyhmmer.easel.Alphabet.amino()
        elif candidate_class in PolynucleotideClass: 
            alphabet = pyhmmer.easel.Alphabet.rna()
        else:
            raise Exception("hmm_produce: Unimplemented candidate class")
                
        seqs_a = muscle_align_N_seq(iter(seed_sequences))
        name   = "{}".format(candidate_class.value)
        HMM    = hmm_create(name, seqs_a, alphabet)

        if not no_cache:
            hmm_cache(HMM)
        return (candidate_class, HMM )

class HMMs():
    # https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.HMM
    """
    STEPS (bottom-up):
    - given a taxonomic id, produce a set of HMMs 
    - each HMM is representative of a polymer class parametrized by a seed MSA 
    - a sequence is searched against all HMMs in the registry
    """


    def __init__(self, tax_id:int, candidate_classes:list[PolymerClass],  no_cache:bool=False, max_seed_seqs:int=5) -> None:

        self.organism_tax_id = tax_id
        self.class_hmms_seed_sequences : dict[PolymerClass, list[SeqRecord]] = {}
        self.class_hmms_registry       : dict[PolymerClass, HMM]             = {}
    
        def __load_seed_sequences(candidate_class:PolymerClass,tax_id:int, max_seed_seqs:int)->Tuple[PolymerClass, list[SeqRecord]]:
            """This is just for housekeeping. Keeping sequences with which the HMMs were seeded as state on the classifier."""
            return (candidate_class, [*fasta_phylogenetic_correction(candidate_class, tax_id, max_n_neighbors=max_seed_seqs)] )

        with concurrent.futures.ThreadPoolExecutor(max_workers=25) as executor:
            loading_futures = []
            loaded_results  = []

            # ! load seqs
            for candidate_class in candidate_classes:
                future = executor.submit(__load_seed_sequences, candidate_class, tax_id, max_seed_seqs)
                loading_futures.append(future)

            for future in concurrent.futures.as_completed(loading_futures): 
                loaded_results.append(future.result())
            # ! populate seqs records
            for (cls, seedseqs) in loaded_results:
                 self.class_hmms_seed_sequences.update({str( cls.value ):seedseqs})
            # pprint("ADDED SEQUENCES TO CLASSIFIER")
            # pprint(self.class_hmms_seed_sequences)

            # ! clean containers
            loaded_results, loading_futures = [],[]
            # ! build hmms
            for candidate_class in candidate_classes:
                future = executor.submit(hmm_produce,candidate_class, tax_id, self.class_hmms_seed_sequences[candidate_class.value], no_cache=no_cache)
                loading_futures.append(future)
            for future in concurrent.futures.as_completed(loading_futures):
                loaded_results.append(future.result())

            # ! populate hmms
            for (cls, hmm) in loaded_results:
                 self.class_hmms_registry.update({cls.value:hmm})


            # pprint(self.class_hmms_registry)
# # !-
#         for candidate_class in candidate_classes:
#             seqs                                       = [*fasta_phylogenetic_correction(candidate_class, tax_id, max_n_neighbors=max_seed_seqs)]
#             self.seed_sequences[candidate_class.value] = seqs
#             self.hmms_registry[candidate_class.value]  = hmm_produce(candidate_class, tax_id, no_cache=no_cache, max_seed_seq=max_seed_seqs)
#             print("Loded HMM: {}".format(candidate_class))
# # !-


    def classify_seq(self, alphabet, target_seq:DigitalSequence)->list[Tuple[PolymerClass, TopHits]]:
        """analogue of `scan`"""
        _ = []
        for (candidate_class, hmm) in self.class_hmms_registry.items():
            result = seq_evaluate_v_hmm(target_seq,alphabet, hmm)
            _.append(( candidate_class,result ))
        return _



    def scan_seq(self, alphabet:pyhmmer.easel.Alphabet, target_seq:SeqRecord):
        """Construct a scan pipeline for the current classifier"""

        # convert seq records to a list of pyhmmer.DigitalSequences
        query_seqs = []
        for seq_record in target_seqs:
            query_seq  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
            query_seqs.append(query_seq.digitize(alphabet))

        scans = list(pyhmmer.hmmscan(query_seqs,[*self.class_hmms_registry.values()] ))
        return scans

    def info(self)->dict:
        """Get info for the constituent HMMs in the current classifier"""
        _ = {
             "seed_seqs": {
                k: [{"seq":str(seqrecord.seq), "tax_id":seqrecord.id} for seqrecord in v] for (k, v) in self.class_hmms_seed_sequences.items()
             },
             "hmms_registry": {
                    k: {"name": v.name.decode(),  "M":v.M, "nseq":v.nseq, "nseq_effective":v.nseq_effective}
                    for (k, v) in self.class_hmms_registry.items()
                },
             }
        return _

class HMMClassifier():

    # organism_scanners : dict[int, HMMs] 
    # chains            : list[Polymer] 
    # alphabet          : pyhmmer.easel.Alphabet
    # candidate_classes : list[PolymerClass]
    bitscore_threshold: int = 35
    # report            : dict

    def __init__(self, chains: list[Polymer], alphabet:pyhmmer.easel.Alphabet, candidate_classes:list[PolymerClass]) -> None:

        self.organism_scanners = {}
        self.chains            = chains
        self.alphabet          = alphabet
        self.report            = {}
        self.candidate_classes = candidate_classes

    def pick_best_hit(self, hits:list[dict]) -> list[PolymerClass]:
        # TODO: look at filtering allunder self.bitscore threshold (while scaling the scores by the [average of the lowest half * 1.5]/[average of the lowest half])
        if len(hits) == 0:
            return []

        sorted_by_biscore = sorted(hits, key=lambda x: x["bitscore"], reverse=True)
        return [ sorted_by_biscore[0]['class_name'] ] if sorted_by_biscore[0]['bitscore'] > self.bitscore_threshold else []


        # return [ sorted(hits, key=lambda x: x["hit.score"])[0]["hit.name"].decode() ] if len(hits) >0 else []

    def classify_chains(self)->None:
        """This is an alternative implementation of `scan_chains` that does not use `hmmscan`. Waiting on https://github.com/althonos/pyhmmer/issues/53 to resolve"""

        for chain in self.chains:
            try:
                organism_taxid = chain.src_organism_ids[0]
                if organism_taxid not in self.organism_scanners:
                    hmmscanner                             = HMMs(organism_taxid, self.candidate_classes, no_cache = True, max_seed_seqs = 5)
                    self.organism_scanners[organism_taxid] = hmmscanner
                else:
                    hmmscanner = self.organism_scanners[organism_taxid]

                self.report[chain.auth_asym_id] = []
                seq_record                      = chain.to_SeqRecord()
                query_seq                       = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq).digitize(self.alphabet)
                for ( candidate_class, tophits ) in hmmscanner.classify_seq(self.alphabet, query_seq):
                    for hit in tophits:
                           d_hit = {
                                "fasta_seed"        : [str(seqrecord.seq) for seqrecord in hmmscanner.class_hmms_seed_sequences[candidate_class]],
                                "seed_organism_ids" : [seqrecord.id for seqrecord in hmmscanner.class_hmms_seed_sequences[candidate_class]],
                                "class_name"        : candidate_class,
                                "evalue"            : hit.evalue,
                                "bitscore"          : hit.score,
                                "target_organism_id": organism_taxid,
                                "consensus"         : hmmscanner.class_hmms_registry[candidate_class].consensus,
                                "domains"           : [( d.score, d.c_evalue, d.env_from, d.env_to ) for d in hit.domains]
                            }
                           self.report[chain.auth_asym_id].append(d_hit)
            except Exception as e:
                print(e)

    def ___scan_chains(self)->None:
        """DEPRECATED: Waiting on https://github.com/althonos/pyhmmer/issues/53 to resolve"""

        for chain in self.chains:
            organism_taxid = chain.src_organism_ids[0]

            # -- pick scanner'|
            if organism_taxid not in self.organism_scanners:
                    hmmscanner                             = HMMs(organism_taxid, self.candidate_classes, no_cache = True, max_seed_seqs = 5)
                    self.organism_scanners[organism_taxid] = hmmscanner
                    pprint(hmmscanner.class_hmms_registry)
            else:
                hmmscanner = self.organism_scanners[organism_taxid]
            # -- pick scanner.|

            self.report[chain.auth_asym_id] = []
            seq_record = chain.to_SeqRecord()
            query_seq  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
            query_seqs = [query_seq.digitize(self.alphabet)]

            # -- convert seq to easel format
            
            # print("Scanning chain {}.{} against {} HMMs".format(chain.parent_rcsb_id, chain.auth_asym_id, len(hmmscanner.class_hmms_registry)))
            # print("seq:", query_seq.sequence)

            # print("Scanning query seqs", [*query_seqs][0].alphabet, [*query_seqs][0].sequence)
            # print("Against hmms", [*hmmscanner.class_hmms_registry.values()])
            
            for scan in list(pyhmmer.hmmscan(query_seqs,[*hmmscanner.class_hmms_registry.values()], self.alphabet, background =pyhmmer.plan7.Background(self.alphabet))):
                for hit in [*scan]:
                   d_hit = {
                    "fasta_seed"        : [str(seqrecord.seq) for seqrecord in hmmscanner.class_hmms_seed_sequences[hit.name.decode()]],
                    "seed_organism_ids" : [seqrecord.id for seqrecord in hmmscanner.class_hmms_seed_sequences[hit.name.decode()]],
                    "class_name"        : hit.name.decode(),                                                             # <- comes from the emitting hmm
                    "evalue"            : hit.evalue,
                    "bitscore"          : hit.score,
                    "target_organism_id": organism_taxid,
                    "consensus"         : hmmscanner.class_hmms_registry[hit.name.decode()].consensus,
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


