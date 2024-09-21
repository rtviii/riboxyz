from pprint import pprint
import sys
from ribctl.lib.schema.types_ribosome import PolymerClass
sys.path.append('/home/rtviii/dev/riboxyz')
from ribctl.lib.libmsa import Fasta
from Bio.SeqRecord import SeqRecord


# Ok let's not complicate things. The prototype is thus:
# take a motif in 7K00 protein, say uL10
# 1. To Record:
# - locus name and its anchor structure/chain/residue
# - project anchor into the MSA of homologs of other structure
# - project anchor into each other structure PAIRWISE
# 2. Backtrack to structural files indices for each structure, record as a row in a table.


s = Fasta.poly_class_all_seq(PolymerClass('uL4'))

type ResidueRange = list[int]


class SequenceProjection_ManyToOne():
    """Given multiple sequences and residue ranges within them, project the ranges onto a single sequence
    """

    source_sequences: dict[str,tuple[SeqRecord,ResidueRange]]
    target_sequence : SeqRecord


    def __init__(self, source_sequences: dict[str,tuple[SeqRecord,ResidueRange]], target_sequence: SeqRecord) -> None:

      








