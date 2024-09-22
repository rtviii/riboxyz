from pprint import pprint
import sys
from pydantic import BaseModel
from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import BiopythonChain
from ribctl.lib.schema.types_binding_site import ResidueSummary
from ribctl.lib.schema.types_ribosome import PolymerClass
sys.path.append('/home/rtviii/dev/riboxyz')
from ribctl.lib.libmsa import Fasta
from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel
from typing import Dict, List


# Ok let's not complicate things. The prototype is thus:
# take a motif in 7K00 protein, say uL10
# 1. To Record:
# - locus name and its anchor structure/chain/residue
# - project anchor into the MSA of homologs of other structure
# - project anchor into each other structure PAIRWISE
# 2. Backtrack to structural files indices for each structure, record as a row in a table.

s = Fasta.poly_class_all_seq(PolymerClass('uL4'))

type ResidiueIndices = list[int]
type RCSB_ID         = str


class SequenceProjection_ManyToOne():
    """Given multiple sequences and residue ranges within them, project the ranges onto a single sequence
    """

    target_chain : BiopythonChain
    target_struct: RCSB_ID
    polymer_class: PolymerClass
    source_chains: dict[RCSB_ID,tuple[BiopythonChain,list[ResidueSummary]]]
    mappings     : dict[RCSB_ID,dict[ResidiueIndices, ResidiueIndices]]
    weights      : dict[int, list[RCSB_ID]]

    def __init__(self,target_seq:BiopythonChain, target_rcsb_id: str,source_seq_pairs:list[tuple[RCSB_ID,BiopythonChain,list[ResidueSummary]]], polymer_class:str) -> None:

        self.target_chain  = target_seq
        self.target_struct = target_rcsb_id
        self.polymer_class = PolymerClass(polymer_class)

        for rcsb_id,chain,indices in source_seq_pairs:
            self.source_chains[rcsb_id] = (chain,indices)

    def project(self):

        # TODO: put this on a threadpool 
        # with concurrent.futures.ProcessPoolExecutor() as executor:
        #     for number, prime in zip(PRIMES, executor.map(is_prime, PRIMES)):
        #         print('%d is prime: %s' % (number, prime))

        for rcsb_id, (source_chain,source_residues) in self.source_chains.items():
            map_motifs(source_chain, self.target_chain, source_residues, self.polymer_class)

        



 


      








