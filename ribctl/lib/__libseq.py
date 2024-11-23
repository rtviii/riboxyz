from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import  Residue
from Bio import pairwise2
from typing import List, Dict, Tuple, Optional
from functools import lru_cache
import re

from ribctl.lib.schema.primitivs import AMINO_ACIDS, NUCLEOTIDES
from ribctl.lib.schema.types_ribosome import ResidueSummary

class SequenceMappingContainer:
    def __init__(self, chain: Chain):
        self.chain = chain
        # Cache the sequence mappings on initialization
        self._init_sequences()
    
    def _init_sequences(self):
        """Initialize all sequence mappings once during construction"""
        self._flat_seq_data = self._compute_flat_sequence()
        self._primary_seq_data = self._compute_primary_sequence()
    
    @lru_cache(maxsize=1)
    def _compute_primary_sequence(self, represent_noncanonical: str = ".") -> Tuple[str, Dict[int, int]]:
        """Compute primary sequence and mapping (cached)"""
        seq = []
        auth_seq_id_to_primary_ix = {}
        
        for ix, residue in enumerate(self.chain.get_residues()):
            resname = residue.resname
            if resname in AMINO_ACIDS:
                seq.append(ResidueSummary.three_letter_code_to_one(resname))
            elif resname in NUCLEOTIDES:
                seq.append(resname)
            else:
                seq.append(represent_noncanonical)
            auth_seq_id_to_primary_ix[residue.get_id()[1]] = ix
            
        return ''.join(seq), auth_seq_id_to_primary_ix

    @lru_cache(maxsize=1)
    def _compute_flat_sequence(self) -> Tuple[str, Dict[int, Residue], Dict[int, int]]:
        """Compute flat sequence and mappings (cached)"""
        seq = []
        flat_index_to_residue_map = {}
        auth_seq_id_to_flat_index_map = {}
        flat_index = 0
        
        for residue in self.chain.get_residues():
            resname = residue.resname
            if resname in AMINO_ACIDS or resname in NUCLEOTIDES:
                seq.append(ResidueSummary.three_letter_code_to_one(resname))
                flat_index_to_residue_map[flat_index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = flat_index
                flat_index += 1
                
        return ''.join(seq), flat_index_to_residue_map, auth_seq_id_to_flat_index_map
    
    @property
    def primary_sequence(self) -> Tuple[str, Dict[int, int]]:
        return self._primary_seq_data
    
    @property
    def flat_sequence(self) -> Tuple[str, Dict[int, Residue], Dict[int, int]]:
        return self._flat_seq_data

class SeqPairwise:
    def __init__(self, sourceseq: str, targetseq: str, source_residues: List[int]):
        self.src     = sourceseq
        self.src_ids = source_residues
        self.tgt     = targetseq
        
        # Compute alignment once
        alignment    = pairwise2.align.globalxx(self.src, self.tgt, one_alignment_only=True)[0]
        self.src_aln = alignment.seqA
        self.tgt_aln = alignment.seqB
        
        # Create position maps first
        self._src_aln_pos_map = self._create_position_map(self.src_aln)
        self._tgt_aln_pos_map = self._create_position_map(self.tgt_aln)
        
        # Then compute indices using the position maps
        self.aligned_ids = self._compute_aligned_indices()
        self.tgt_ids     = self._compute_target_indices()
    
    def _create_position_map(self, seq: str) -> Dict[int, int]:
        """Create a mapping of position to original index (excluding gaps)"""
        pos_map = {}
        orig_pos = 0
        for i, char in enumerate(seq):
            if char != '-':
                pos_map[i] = orig_pos
                orig_pos += 1
        return pos_map
    
    def _compute_aligned_indices(self) -> List[int]:
        """Compute all aligned indices in one pass"""
        indices = []
        for src_resid in self.src_ids:
            orig_count = 0
            for i, char in enumerate(self.src_aln):
                if char != '-':
                    if orig_count == src_resid:
                        indices.append(i)
                        break
                    orig_count += 1
        return indices
    
    def _compute_target_indices(self) -> List[int]:
        """Compute all target indices in one pass"""
        indices = []
        for aln_resid in self.aligned_ids:
            if self.tgt_aln[aln_resid] != '-':
                orig_pos = self._tgt_aln_pos_map.get(aln_resid)
                if orig_pos is not None:
                    indices.append(orig_pos)
        return indices

    @staticmethod
    def highlight_indices(sequence: str, ixs: List[int], color: int = 91) -> str:
        """Optimized highlighting using a single join operation"""
        highlighted = []
        for i, v in enumerate(sequence):
            if i in ixs:
                highlighted.append(f"\033[{color}m{v}\033[0m")
            else:
                highlighted.append(v)
        return ''.join(highlighted)


def map_motifs(
    source_chain: SequenceMappingContainer,
    target_chain: SequenceMappingContainer,
    bound_residues: List[ResidueSummary],
    polymer_class: str,
    verbose: bool = False
) -> Tuple[str, str, List[Residue]]:
    
    # Get cached sequences and mappings
    src_flat_seq, src_flat_to_res, src_auth_to_flat = source_chain.flat_sequence
    tgt_flat_seq, tgt_flat_to_res, _                = target_chain.flat_sequence
    
    # Filter bound residues once
    valid_residues = {*NUCLEOTIDES, *AMINO_ACIDS.keys()}
    src_bound_flat_indices = [ src_auth_to_flat[res.auth_seq_id] for res in bound_residues if res.label_comp_id in valid_residues ]
    
    # Perform sequence alignment and mapping
    mapper             = SeqPairwise(src_flat_seq, tgt_flat_seq, src_bound_flat_indices)
    tgt_bound_residues = [tgt_flat_to_res[idx] for idx in mapper.tgt_ids]
    
    if verbose:
        primary_seq_source, auth_to_primary_source = source_chain.primary_sequence
        primary_seq_target, auth_to_primary_target = target_chain.primary_sequence
        
        print(f"\n\n\t\t [{polymer_class}] ")
        print("[\033[95mSource\033[0m Primary]\t", 
              SeqPairwise.highlight_indices(
                  primary_seq_source,
                  [auth_to_primary_source[res.auth_seq_id] for res in bound_residues]
              ))
        print("[\033[95mSource\033[0m Flat   ]\t",
              SeqPairwise.highlight_indices(src_flat_seq, src_bound_flat_indices))
        print("[\033[95mSource\033[0m Aligned]\t",
              mapper.highlight_indices(mapper.src_aln, mapper.aligned_ids))
        print("[\033[96mTarget\033[0m Aligned]\t",
              mapper.highlight_indices(mapper.tgt_aln, mapper.aligned_ids))
        print("[\033[96mTarget\033[0m Flat   ]\t",
              SeqPairwise.highlight_indices(tgt_flat_seq, mapper.tgt_ids))
        print("[\033[96mTarget\033[0m Primary]\t",
              SeqPairwise.highlight_indices(
                  primary_seq_target,
                  [auth_to_primary_target[res.get_id()[1]] for res in tgt_bound_residues]
              ))
    
    return source_chain.primary_sequence[0], target_chain.primary_sequence[0], tgt_bound_residues