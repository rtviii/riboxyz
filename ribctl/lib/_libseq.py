import sys

from ribctl.lib.schema.primitivs import AMINO_ACIDS, NUCLEOTIDES
from ribctl.lib.schema.types_ribosome import ResidueSummary
sys.path.append("/home/rtviii/dev/riboxyz")
from typing import Optional
import re
from typing import List
import warnings
from Bio import ( BiopythonDeprecationWarning, )
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import pairwise2
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain

# lets say you align 4ug0.A to 4u3m.F
# 1. take biopython struct of 4ugo, get chain A
# 2. take biopython struct of 4u3m, get chain F
# 3. wrap both in SequenceMappingContainer
# 4. use flat sequence for alignment

class SequenceMappingContainer(Chain):
    """ 
    A container for keeping track of the correspondence between the structural and sequence indices within a give polymer chain.
    \nStructural data(`mmcif`) frequently has unresolved and modified residues adheres to the author's numbering (which is arbitrary for all intents and purposes, ex. starts at, say, 7).
    More here:https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
    
    This container keeps pointers between the naive (convenient) indices and the structural residues. 
    I refer to the following in the code:

    - **auth_seq_id** is the author-assigned residue number and is frequently used to refer to structural components. (ex. in Molstar)
    - **flat_index** is the the straightforward arithmetic (0-start) numbering. I produce it by removing all the modified residues from the "primary sequence"
    - **primary_sequence** is the sequence of structural residues (as packed into the BioPython `Chain` object) represented as a string.
    - **flat_sequence** is the **primary_sequence** with the modified residues removed, represented as a string.
    
    There is lots to optimize in this code (it builds index->Residue<object> maps by enumeration),
    but ideally this is taken care of at the parser level or at the deposition level.

    Again, see more: 

    - https://proteopedia.org/wiki/index.php/Unusual_sequence_numbering
    - https://bioinformatics.stackexchange.com/questions/14210/pdb-residue-numbering
    - https://bioinformatics.stackexchange.com/questions/20458/how-is-the-canonical-version-entity-poly-pdbx-seq-one-letter-code-obtaine
    - https://www.biostars.org/p/9588718/
    - https://stackoverflow.com/questions/45466408/biopython-resseq-doesnt-match-pdb-file
    
    """

    chain                        : Chain

    # from flat to residue that contributes the index to primary
    flat_index_to_residue_map    : dict[int, Residue]

    # from primary to flat
    auth_seq_id_to_flat_index_map: dict[int, int]

    @property
    def primary_sequence(self, represent_noncanonical_as:str=".") -> tuple[str, dict[int,int]]:
        seq = ""
        auth_seq_id_to_primary_ix = {}
        for ix, residue in enumerate(self.chain.get_residues()):
            if residue.resname in [*AMINO_ACIDS.keys()]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
            elif residue.resname in [*NUCLEOTIDES]:
                seq = seq + residue.resname
            else:
                seq = seq + represent_noncanonical_as
            auth_seq_id_to_primary_ix[residue.get_id()[1]] = ix

        return seq, auth_seq_id_to_primary_ix

    @property
    def flat_sequence(self) -> tuple[str, dict[int,Residue], dict[int,int]]:
        res: list[Residue] = [*self.chain.get_residues()]

        flat_index_to_residue_map     = {}
        auth_seq_id_to_flat_index_map = {}
        seq                           = ""
        flat_index                    = 0

        for residue in res:
            if residue.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                flat_index_to_residue_map[flat_index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = flat_index
                flat_index += 1
            else:
                continue
        return seq, flat_index_to_residue_map, auth_seq_id_to_flat_index_map

    def __init__(self, chain: Chain):
        self.chain = chain

class SeqPairwise:
    def __init__(self, sourceseq: str, targetseq: str, source_residues: list[int]):
        """A container for origin and target sequences when matching residue indices in the source sequence to the target sequence.
         - return the map {int:int} between the two sequences
                 CASE 1. The first one is longer:
                     Initial

        Common ix    0 1 2 3 4 5 6 7 8 9             Aligned ix   0 1 2 3 4 5 6 7 8 9
        Canonical    X Y G G H A S D S D    ----->   Canonical    X Y G G H A S D S D
        Structure    Y G G H A S D                   Structure    - Y G G H A S D - -

         IMPORTANT: BOTH SEQUENCES ARE ASSUMED TO HAVE NO GAPS ( at least not represeneted as "-"). That will screw up the arithmetic.
        """

        # *  indices of the given residues in the source sequence.
        self.src    : str       = sourceseq
        self.src_ids: list[int] = source_residues

        # * Indices of the corresponding residues in target sequence. To be filled.
        self.tgt: str = targetseq
        self.tgt_ids: list[int] = []

        _ = pairwise2.align.globalxx(self.src, self.tgt, one_alignment_only=True)

        self.src_aln = _[0].seqA
        self.tgt_aln = _[0].seqB

        self.aligned_ids = []

        for src_resid in self.src_ids:
            self.aligned_ids.append(self.forwards_match(self.src_aln, src_resid))

        for aln_resid in self.aligned_ids:
            tgt_aln_index = self.backwards_match(self.tgt_aln, aln_resid)
            if tgt_aln_index == None:
                continue
            else:
                self.tgt_ids.append(tgt_aln_index)

    def forwards_match(
        self, aligned_source_sequence: str, original_residue_index: int
    ) -> int:
        """Returns the index of a source-sequence residue in the aligned source sequence. Basically, "count forward including gaps" """
        if original_residue_index > len(aligned_source_sequence):
            raise IndexError( f"Passed residue with invalid index ({original_residue_index}) to back-match to target.Seqlen aligned:{len(aligned_source_sequence)}" )

        original_residues_count = 0
        for aligned_ix, char in enumerate(aligned_source_sequence):
            if original_residues_count == original_residue_index:
                if char == "-":
                    continue
                else:
                    return aligned_ix
            if char == "-":
                continue
            else:
                original_residues_count += 1

        raise ValueError(
            f"Residue with index {original_residue_index} not found in the aligned source sequence after full search. Logical errory, likely."
        )

    def backwards_match(
        self, aligned_target_sequence: str, aligned_residue_index: int
    ) -> int | None:
        """Returns the target-sequence index of a residue in the [aligned] target sequence. Basically, "count back ignoring gaps" """
        if aligned_residue_index > len(aligned_target_sequence):
            raise IndexError( f"Passed residue with invalid index ({aligned_residue_index}) to back-match to target.Seqlen:{len(aligned_target_sequence)}" )

        if aligned_target_sequence[aligned_residue_index] == "-":
            return None

        original_residues_index = 0
        for aligned_ix, char in enumerate(aligned_target_sequence):
            if aligned_ix == aligned_residue_index:
                return original_residues_index
            if char == "-":
                continue
            else:
                original_residues_index += 1

    @staticmethod
    def highlight_subseq(sequence: str, subsequence: str, index: int = None):
        """Highlight subsequence"""
        CRED = "\033[91m"
        CEND = "\033[0m"
        _ = []
        if index != None:
            return (
                sequence[: index - 1]
                + CRED
                + sequence[index]
                + CEND
                + sequence[index + 1 :]
            )
        for item in re.split(re.compile(f"({subsequence})"), sequence):
            if item == subsequence:
                _.append(CRED + item + CEND)
            else:
                _.append(item)
        return "".join(_)

    @staticmethod
    def highlight_indices(sequence: str, ixs: List[int], color: int = 91):
        """Highlight indices"""
        CRED = "\033[{}m".format(color)
        CEND = "\033[0m"
        _ = ""
        for i, v in enumerate(sequence):
            if i in ixs:
                _ += CRED + v + CEND
            else:
                _ += v
        return _


def map_motifs(source_chain:SequenceMappingContainer, target_chain:SequenceMappingContainer, bound_residues:list[ResidueSummary], polymer_class:str, verbose:bool=False)->tuple[str,str,list[Residue]]:

    bpchain_source = source_chain
    bpchain_target = target_chain

    #! SOURCE & TARGET MAPS
    [ src_flat_structural_seq, src_flat_idx_to_residue_map, src_auth_seq_id_to_flat_index_map, ] = bpchain_source.flat_sequence
    [ tgt_flat_structural_seq, tgt_flat_idx_to_residue_map, tgt_auth_seq_id_to_flat_index_map, ] = bpchain_target.flat_sequence

    # ! Bound residues  [in STRUCTURE SPACE]
    src_bound_auth_seq_idx = [ (residue.auth_seq_id, residue.label_comp_id) for residue in bound_residues ]

    primary_seq_source, auth_seq_to_primary_ix_source = ( bpchain_source.primary_sequence )
    primary_seq_target, auth_seq_to_primary_ix_target = ( bpchain_target.primary_sequence )


    src_bound_flat_indices = [ src_auth_seq_id_to_flat_index_map[index] for index, label in filter( lambda x: x[1] in [*NUCLEOTIDES, *AMINO_ACIDS.keys()], src_bound_auth_seq_idx, ) ]

    M = SeqPairwise( src_flat_structural_seq, tgt_flat_structural_seq, src_bound_flat_indices )

    tgt_bound_flat_indices = M.tgt_ids
    tgt_bound_residues = [ tgt_flat_idx_to_residue_map[idx] for idx in tgt_bound_flat_indices ]

    if verbose:
        print("\n\n\t\t [{}] ".format(polymer_class))
        print( "[\033[95mSource\033[0m Primary]\t", SeqPairwise.highlight_indices( primary_seq_source, [ auth_seq_to_primary_ix_source[index] for index, label in src_bound_auth_seq_idx ]))
        print( "[\033[95mSource\033[0m Flat   ]\t", SeqPairwise.highlight_indices(src_flat_structural_seq, src_bound_flat_indices), )
        print( "[\033[95mSource\033[0m Aligned]\t", M.highlight_indices(M.src_aln, ixs=M.aligned_ids) )
        print( "[\033[96mTarget\033[0m Aligned]\t", M.highlight_indices(M.tgt_aln, ixs=M.aligned_ids) )
        print( "[\033[96mTarget\033[0m Flat   ]\t", SeqPairwise.highlight_indices(tgt_flat_structural_seq, tgt_bound_flat_indices), )
        print( "[\033[96mTarget\033[0m Primary]\t", SeqPairwise.highlight_indices( primary_seq_target, [ auth_seq_to_primary_ix_target[residue.get_id()[1]] for residue in tgt_bound_residues ], ), )

    return primary_seq_source, primary_seq_target, tgt_bound_residues