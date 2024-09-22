import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from typing import Optional
import re
from typing import List
import warnings
from Bio import ( BiopythonDeprecationWarning, )
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import pairwise2
from ribctl.lib.schema.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
)

from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.Seq import Seq
from ribctl.lib.schema.types_binding_site import ( ResidueSummary, )

#! Transposition methods
class BiopythonChain(Chain):

    chain                        : Chain
    flat_index_to_residue_map    : dict[int, Residue]
    auth_seq_id_to_flat_index_map: dict[int, int]

    def __init__(self, chain: Chain):
        self.chain = chain

    @property
    def primary_sequence(self) -> tuple[str, dict]:

        represent_noncanonical_as: Optional[str] = "."
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
    def flat_sequence(self) -> tuple[str, dict, dict]:
        res: list[Residue] = [*self.chain.get_residues()]

        flat_index_to_residue_map = {}
        auth_seq_id_to_flat_index_map = {}
        seq = ""
        flat_index = 0
        for residue in res:
            if residue.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                flat_index_to_residue_map[flat_index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = flat_index
                flat_index += 1
            else:
                continue
        return seq, flat_index_to_residue_map, auth_seq_id_to_flat_index_map

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
            raise IndexError(
                f"Passed residue with invalid index ({original_residue_index}) to back-match to target.Seqlen aligned:{len(aligned_source_sequence)}"
            )
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
            raise IndexError(
                f"Passed residue with invalid index ({aligned_residue_index}) to back-match to target.Seqlen:{len(aligned_target_sequence)}"
            )

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
    def hl_subseq(sequence: str, subsequence: str, index: int = None):
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
    def hl_ixs(sequence: str, ixs: List[int], color: int = 91):
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

class SeqMap:

    mapping: dict[int, int]

    seq_canonical: str
    seq_structural: str

    seq_canonical_aligned: str
    seq_structural_aligned: str

    def __init__(self, canonical: str, structure: str):
        self.seq_canonical  = canonical
        self.seq_structural = structure

        alignments = pairwise2.align.globalxx(Seq(canonical), Seq(structure))
        aligned_canonical, aligned_structure = alignments[0][0], alignments[0][1]

        self.seq_canonical_aligned  = aligned_canonical
        self.seq_structural_aligned = aligned_structure

        mapping = {}

        # print("inspecting")
        # print(self.seq_canonical_aligned)
        # print(self.seq_structural_aligned)
        # print(*zip(aligned_canonical, aligned_structure))

        canonical_index = 0
        structure_index = 0
        for canonical_char, structural_char in zip(
            aligned_canonical, aligned_structure
        ):
            if canonical_char != "-":
                if structural_char != "-":
                    mapping[canonical_index] = structure_index
                    structure_index += 1
                else:
                    mapping[canonical_index] = -1
                canonical_index += 1
            elif canonical_char == "-":
                continue
                # warnings.warn(f"Unexpected gap in canonical sequence at aligned position {canonical_index}. This shouldn't happen with the original canonical sequence.")

        self.mapping = mapping

    def retrieve_index(self, key: int) -> int | None:
        "Get the STRUCTURAL sequence index corresponding to the CANONICAL sequence index <key> if any, otherwise None"
        if key not in self.mapping:
            raise KeyError(f"Key {key} not found in mapping")
        if self.mapping[key] == -1:
            return None
        return self.mapping[key]

    def retrieve_motif(self, keys: list[int]) -> tuple[str, str]:
        can_subseq    = ""
        struct_subseq = ""

        for i in keys:
            can_subseq = can_subseq + self.seq_canonical[i]
            struct_index = self.retrieve_index(i)
            if struct_index == None:
                struct_subseq = struct_subseq + "-"
            elif struct_index != None:
                struct_subseq = struct_subseq + self.seq_structural[struct_index]
        if struct_subseq == "" or list(set(list(struct_subseq)))[0] == "-":
            raise ValueError( "No structural sequence found for the given canonical sequence" )
        return can_subseq, struct_subseq