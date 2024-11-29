# Create Fasta object with sequences
import os
from pprint import pprint
from ribctl import ASSETS
from ribctl.lib.libmsa import Fasta, FastaBalancer
from ribctl.lib.libseq import SequenceMappingContainer, get_conservation_scores, map_alignment_to_structure
from ribctl.ribosome_ops import RibosomeOps


fb               = FastaBalancer(Fasta(os.path.join(ASSETS['fasta_proteins_cytosolic'], 'uL4.fasta')))
fasta_obj        = fb.balance_dataset_sparse(30)
alignment_scores = get_conservation_scores(fasta_obj)
# struct           = RibosomeOps('4UG0').assets.biopython_structure()[0]
# polymer          = RibosomeOps('4UG0').get_poly_by_polyclass('uL4')
chain = RibosomeOps('4UG0').get_biopython_chain_by_polymer_class('uL4')


structure_scores = map_alignment_to_structure(alignment_scores,SequenceMappingContainer(chain))
pprint(structure_scores)