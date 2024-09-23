from pprint import pprint
import sys
from pydantic import BaseModel
from Bio.PDB.Residue import Residue
from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import BiopythonChain
from ribctl.lib.schema.types_binding_site import ResidueSummary
from ribctl.lib.schema.types_ribosome import PolymerClass
from functools import partial
sys.path.append('/home/rtviii/dev/riboxyz')
from ribctl.lib.libmsa import Fasta
from Bio.SeqRecord import SeqRecord
from pydantic import BaseModel
from typing import Dict, List
from concurrent.futures import ProcessPoolExecutor


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

    target_chain   : BiopythonChain
    target_struct  : RCSB_ID
    polymer_class  : PolymerClass

    source_chains  : dict[RCSB_ID,tuple[BiopythonChain,list[ResidueSummary]]]
    target_mappings: dict[RCSB_ID,list[Residue]]

    weights        : dict[int, list[RCSB_ID]]

    def __init__(self,target_seq:BiopythonChain, target_rcsb_id: str,source_seq_pairs:list[tuple[RCSB_ID,BiopythonChain,list[ResidueSummary]]], polymer_class:str) -> None:

        self.target_chain  = target_seq
        self.target_struct = target_rcsb_id
        self.polymer_class = PolymerClass(polymer_class)

        for rcsb_id,chain,indices in source_seq_pairs:
            self.source_chains[rcsb_id] = (chain,indices)

    def project(self):
        def map_motifs_wrapper(args):
            source_chain, target_chain, source_residues, polymer_class = args
            return map_motifs(source_chain, target_chain, source_residues, polymer_class)

        futures = []
        with ProcessPoolExecutor(max_workers=10) as executor:
            for rcsb_id, (source_chain, source_residues) in self.source_chains.items():
                    partial_map_motifs = partial(map_motifs_wrapper, (source_chain, self.target_chain, source_residues, self.polymer_class))
                    future             = executor.submit(partial_map_motifs)
                    futures.append((rcsb_id, future))


        for rcsb_id, future in futures:
            _s, _t, tgt_bound_residues = future.result()
            self.target_mappings[rcsb_id] = tgt_bound_residues


        # for rcsb_id, (source_chain,source_residues) in self.source_chains.items():
        #     _s, _t, tgt_bound_residues    = map_motifs(source_chain, self.target_chain, source_residues, self.polymer_class)
        #     self.target_mappings[rcsb_id] = tgt_bound_residues

    def get_weights(self):
        self.target_mappings
        ...


