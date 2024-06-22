import json
import os
import pickle
from pprint import pprint
from typing import  Iterator,  Tuple
from Bio.SeqRecord import SeqRecord
import pyhmmer
from ribctl import ASSETS
from ribctl.lib.libmsa import Fasta, muscle_align_N_seq, phylogenetic_neighborhood
from ribctl.lib.schema.types_ribosome import RNA, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, PolynucleotideClass, PolypeptideClass, Protein, CytosolicProteinClass, CytosolicProteinClass, tRNA
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock, DigitalSequence
from pyhmmer.plan7 import Pipeline, HMM , TopHits
from ribctl.logs.loggers import get_classification_logger

import concurrent.futures

# logger= get_classification_logger()
hmm_cachedir = ASSETS['__hmm_cache']





def digitize_seq_record(seq_record:SeqRecord, alphabet:Alphabet)->DigitalSequence:
    """Convert a SeqRecord to a DigitalSequence"""
    seq_  = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=seq_record.seq)
    return seq_.digitize(alphabet)


def pick_best_hmm_hit(matches_dict:dict[PolymerClass, list[float]], chain_info:Polymer)->PolymerClass | None:
    """Given a dictionary of sequence-HMMe e-values, pick the best candidate class"""
    results = []
    # if len([ item for x in list(matches_dict.values()) for item in x ]) > 0:
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
            logger.warning("{}.{} : Multiple sensible matches detected, picked {} (smallest e-value) : \n {}".format(chain_info.parent_rcsb_id, chain_info.auth_asym_id,results[0][0], results))

    return results[0]['candidate_class']


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

    prot_classes = [*list(CytosolicProteinClass),*list(LifecycleFactorClass),*list(MitochondrialProteinClass )]
    rna_classes  = [*list(PolynucleotideClass)]

    if candidate_class in prot_classes:
        alphabet = pyhmmer.easel.Alphabet.amino()
    elif candidate_class in rna_classes: 
        alphabet = pyhmmer.easel.Alphabet.rna()
    else:
        raise Exception("hmm_produce: Unimplemented candidate class")
            
    seqs_a = muscle_align_N_seq(iter(seed_sequences))
    name   = "{}".format(candidate_class.value)
    HMM    = hmm_create(name, seqs_a, alphabet)

    return (candidate_class, HMM )

def _obtain_phylogenetic_nbhd_task(base_taxids:list[int], fasta_record:Fasta, polymer_class:PolymerClass, taxid:int, max_n_neighbors:int):
    phylo_nbhd_ids = phylogenetic_neighborhood(base_taxids, str(taxid), max_n_neighbors)
    seqs           = fasta_record.pick_taxids(phylo_nbhd_ids)
    return polymer_class, iter(seqs)

class PolymerClassFastaRegistry():
    """This is a wrapper around fasta records for all polymer classes. Only a few sequences should be picked from it at a time."""
    registry_fasta       : dict[PolymerClass, Fasta ]     = {}
    registry_all_tax_ids : dict[PolymerClass, list[int] ] = {}

    def __init__(self) :

        for candidate_class in CytosolicProteinClass:
            fasta_path                             = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")
            self.registry_fasta[candidate_class]   = Fasta(fasta_path)
            self.registry_all_tax_ids[candidate_class] = Fasta(fasta_path).all_taxids("int")

        for candidate_class in MitochondrialProteinClass:
            fasta_path                             = os.path.join(ASSETS["fasta_proteins_mitochondrial"], f"{candidate_class.value}.fasta")
            self.registry_fasta[candidate_class]   = Fasta(fasta_path)
            self.registry_all_tax_ids[candidate_class] = Fasta(fasta_path).all_taxids("int")

        for candidate_class in PolynucleotideClass:
            fasta_path                             = os.path.join(ASSETS["fasta_rna"], f"{candidate_class.value}.fasta")
            self.registry_fasta[candidate_class]   = Fasta(fasta_path)
            self.registry_all_tax_ids[candidate_class] = Fasta(fasta_path).all_taxids("int")

        for candidate_class in ElongationFactorClass:
            fasta_path                             = os.path.join(ASSETS["fasta_factors_elongation"], f"{candidate_class.value}.fasta")
            self.registry_fasta[candidate_class]   = Fasta(fasta_path)
            self.registry_all_tax_ids[candidate_class] = Fasta(fasta_path).all_taxids("int")

        for candidate_class in InitiationFactorClass:
            fasta_path                             = os.path.join(ASSETS["fasta_factors_initiation"], f"{candidate_class.value}.fasta")
            self.registry_fasta[candidate_class]   = Fasta(fasta_path)
            self.registry_all_tax_ids[candidate_class] = Fasta(fasta_path).all_taxids("int")
            
        #-------------------- Initated all polymer class sequences. Check that they are not empty.
        
        for (polymer_class, fasta_record) in self.registry_fasta.items():
            if len( list(fasta_record.records) ) == 0 :
                raise Exception("Empty fasta file: {}".format(fasta_record))

    def get_seed_sequences_for_taxid(self,taxid:int, max_n_neighbors:int = 5)->dict[PolymerClass, list[SeqRecord]]:
        _ = {}
        with concurrent.futures.ProcessPoolExecutor(max_workers=15) as executor:
            futures = []
            for (polymer_class, fasta_record) in self.registry_fasta.items():
                futures.append(executor.submit(
                                    _obtain_phylogenetic_nbhd_task,
                                    self.registry_all_tax_ids[polymer_class],
                                    fasta_record,
                                    polymer_class,
                                    taxid,
                                    max_n_neighbors))
        for completed_future in concurrent.futures.as_completed(futures):
                    polymer_class, seqs  = completed_future.result()
                    _[polymer_class] = seqs
        return _ 

class PolymerClassesOrganismScanner():
    # https://pyhmmer.readthedocs.io/en/stable/api/plan7.html#pyhmmer.plan7.HMM
    """
    STEPS (bottom-up):
    - given a taxonomic id, produce a set of HMMs 
    - each HMM is representative of a polymer class parametrized by a seed MSA 
    - a sequence is searched against all HMMs in the registry
    """

    filename_hmms_registry:str
    filename_seed_seqs:str


    @property
    def hmms_amino(self)->list[tuple[PolymerClass,HMM]]:
        return list(filter(lambda x: x[0] in [*list(CytosolicProteinClass),*list(LifecycleFactorClass),*list(MitochondrialProteinClass )], self.hmms_registry.items()))

    @property
    def hmms_rna(self)->list[tuple[PolymerClass,HMM]]:
        return list(filter(lambda x: x[0] in [*list(PolynucleotideClass)], self.hmms_registry.items()))

    def cache_hmms(self):
        with open(os.path.join(self.filename_hmms_registry), 'wb')  as outfile:
            pickle.dump(self.hmms_registry, outfile)
            print("Saved {}".format(self.filename_hmms_registry))

    def cache_seed_seqs(self):
        with open(os.path.join(self.filename_seed_seqs), 'wb')  as outfile:
            pickle.dump(self.seed_sequences, outfile)
            print("Saved {}".format(self.filename_hmms_registry))

    def cache_load_hmms(self):
        if os.path.exists(self.filename_hmms_registry):
            with open(self.filename_hmms_registry, 'rb')  as infile:
                self.hmms_registry = pickle.load(infile)
                print("Opened {}".format(self.filename_hmms_registry))
                # pprint(self.hmms_registry)
           
    def cache_load_seed_seqs(self):
        if os.path.exists(self.filename_seed_seqs):
            with open(self.filename_seed_seqs, 'rb')  as infile:
                self.seed_sequences = pickle.load(infile)
                print("Opened {}".format(self.filename_seed_seqs))
                # pprint(self.seed_sequences)

    def __init__(self, organism_taxid:int, no_cache:bool=False, max_seed_seqs:int=5) -> None:

        print("Initializing PolymerClassesOrganismScanner for taxid: {}".format(organism_taxid))

        self.organism_tax_id:int = organism_taxid
        self.filename_hmms_registry       = os.path.join(ASSETS['__hmm_cache'], "taxid_{}_hmm_scanner".format(organism_taxid))
        self.filename_seed_seqs = os.path.join(ASSETS['__hmm_cache'], "taxid_{}_seed_seqs".format(organism_taxid))

        self.seed_sequences: dict[PolymerClass, list[SeqRecord]] = {}
        self.hmms_registry : dict[PolymerClass, HMM]             = {}

        self.cache_load_hmms()
        self.cache_load_seed_seqs()
        

        if ( self.hmms_registry != {} and self.seed_sequences != {} ):
            return

        self.seed_sequences =  PolymerClassFastaRegistry().get_seed_sequences_for_taxid(organism_taxid, max_n_neighbors=max_seed_seqs)

        futures = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
            for ( candidate_class, seqs ) in self.seed_sequences.items():
                future = executor.submit(hmm_produce,candidate_class, organism_taxid, seqs, no_cache=no_cache)
                futures.append(future)

        for completed_future in concurrent.futures.as_completed(futures):
             poly_class, hmm = completed_future.result()
             self.hmms_registry.update({poly_class.value:hmm})
        
        if not no_cache:
            self.cache_hmms()
            self.cache_seed_seqs()
            

    def classify_seq(self, alphabet, target_seq:DigitalSequence)->list[Tuple[PolymerClass, TopHits]]:

        def seq_evaluate_v_hmm(seq:DigitalSequence,alphabet:Alphabet, hmm:HMM, T:int=35)->TopHits:
            """Fit a sequence to a given HMM"""
            dsb   = DigitalSequenceBlock(alphabet, [seq])
            return pyhmmer.plan7.Pipeline(alphabet=alphabet, T=T).search_hmm(hmm,dsb)
        _ = []

        hmms = self.hmms_amino if alphabet == pyhmmer.easel.Alphabet.amino() else self.hmms_rna
        for (candidate_class, hmm) in hmms:
            result = seq_evaluate_v_hmm(target_seq,alphabet, hmm)
            _.append(( candidate_class,result ))
        return _

    def info(self)->dict:
        """Get info for the constituent HMMs in the current classifier"""
        _ = {
             "seed_seqs": {
                k: [{"seq":str(seqrecord.seq), "tax_id":seqrecord.id} for seqrecord in v] for (k, v) in self.seed_sequences.items()
             },
             "hmms_registry": {
                    k: {"name": v.name.decode(),  "M":v.M, "nseq":v.nseq, "nseq_effective":v.nseq_effective}
                    for (k, v) in self.hmms_registry.items()
                },
             }
        return _

class HMMClassifier():
    """
    HMMClassifier takes a set of polymer chains and "candidate" polymer classes and returns a classification report for which polymer class best matches the given polymer chain.

    In a nutshell it's an elaborate search engine: 
    - each polymer is scanned against all classes' HMM's.

    Each class'es HMM is seeded with a set of known sequences of representatives of that class of [ 5 ] "phylogenetically close" organisms.

    Given that among the list of target polymers multiple taxonomic ids might be present (unusual), we keep the the "organism_scanners" field.
    If a new taxonomic id appears, a new "scanner" is created and stored

    
    """

    organism_scanners : dict[int, PolymerClassesOrganismScanner]
    chains            : list[Polymer]
    alphabet          : pyhmmer.easel.Alphabet
    candidate_classes : list[PolymerClass]
    bitscore_threshold: int = 35
    report            : dict

    def __init__(self, chains: list[Polymer], alphabet:pyhmmer.easel.Alphabet, candidate_classes:list[PolymerClass]) -> None:
        self.alphabet          = alphabet
        self.candidate_classes = candidate_classes
        self.chains            = chains
        self.organism_scanners = {}
        self.report            = {}

    def pick_best_hit(self, hits:list[dict]) -> list[PolymerClass]:
        if len(hits) == 0:
            return []

        sorted_by_biscore = sorted(hits, key=lambda x: x["bitscore"], reverse=True)
        return [ sorted_by_biscore[0]['class_name'] ] if sorted_by_biscore[0]['bitscore'] > self.bitscore_threshold else []

    def classify_chains(self)->None:
        for chain in self.chains:
                # --- If the the scanner for this taxid is present, use it. Otherwise create it.
                if len(chain.src_organism_ids)< 1:
                    continue
                organism_taxid = chain.src_organism_ids[0]
                if organism_taxid not in self.organism_scanners:
                    hmmscanner                             = PolymerClassesOrganismScanner(organism_taxid, max_seed_seqs = 5)
                    self.organism_scanners[organism_taxid] = hmmscanner
                else:
                    hmmscanner = self.organism_scanners[organism_taxid]
                # --- If the the scanner for this taxid is present, use it. Otherwise create it.

                self.report[chain.auth_asym_id] = []
                seq_record                      = chain.to_SeqRecord()
                query_seq                       = pyhmmer.easel.TextSequence(name=bytes(seq_record.id,'utf-8'), sequence=str( seq_record.seq )).digitize(self.alphabet)
                cls_hits_tuples                 = hmmscanner.classify_seq(self.alphabet, query_seq)

                for ( candidate_class, tophits ) in cls_hits_tuples:
                    for hit in tophits:
                           d_hit = {
                                "fasta_seed"        : [str(seqrecord.seq) for seqrecord in hmmscanner.seed_sequences[candidate_class]],
                                "seed_organism_ids" : [seqrecord.id for seqrecord in hmmscanner.seed_sequences[candidate_class]],
                                "class_name"        : candidate_class,
                                "evalue"            : hit.evalue,
                                "bitscore"          : hit.score,
                                "target_organism_id": organism_taxid,
                                "consensus"         : hmmscanner.hmms_registry[candidate_class].consensus,
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


