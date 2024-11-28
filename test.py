# Create Fasta object with sequences
import os
from pprint import pprint
from ribctl import ASSETS
from ribctl.lib.libmsa import Fasta, FastaBalancer
from ribctl.lib.libseq import SequenceMappingContainer, get_conservation_scores, map_alignment_to_structure
from ribctl.ribosome_ops import RibosomeOps


fb= FastaBalancer(Fasta(os.path.join(ASSETS['fasta_proteins_cytosolic'], 'uL4.fasta')))

fasta_obj = fb.balance_dataset_sparse(30)

# Get alignment conservation scores
alignment_scores = get_conservation_scores(fasta_obj)

# Map to structure (if needed)
struct = RibosomeOps('4UG0').assets.biopython_structure()[0]
polymer = RibosomeOps('4UG0').get_poly_by_polyclass('uL4')
if polymer == None:
    exit("Polymer not found")

for chain in struct:
    if chain.id == polymer.auth_asym_id:
        structure_scores = map_alignment_to_structure(alignment_scores,         SequenceMappingContainer(chain))
        pprint(structure_scores)
        break