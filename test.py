# Create Fasta object with sequences
import json
import os
from pprint import pprint
from ribctl import ASSETS, ASSETS_PATH
from ribctl.lib.libmsa import Fasta, FastaBalancer
from ribctl.lib.libseq import (
    SequenceMappingContainer,
    get_conservation_scores,
    map_alignment_to_structure,
)
from ribctl.ribosome_ops import RibosomeOps


# # fb               = FastaBalancer(Fasta(os.path.join(ASSETS['fasta_proteins_cytosolic'], 'uL4.fasta')))
# # fasta_obj        = fb.balance_dataset_sparse(30)
# # alignment_scores = get_conservation_scores(fasta_obj)
# # struct           = RibosomeOps('4UG0').assets.biopython_structure()[0]
# # polymer          = RibosomeOps('4UG0').get_poly_by_polyclass('uL4')
# chain = RibosomeOps('7K00').get_biopython_chain_by_polymer_class('23SrRNA')

# seq, mpa= SequenceMappingContainer(chain).primary_sequence

# print(seq.rstrip('.'))
# print(len(seq.rstrip('.')))

# structure_scores = map_alignment_to_structure(alignment_scores,SequenceMappingContainer(chain))
# pprint(structure_scores)


file = os.path.join(ASSETS_PATH, "ligands", "ligand_classes.json")
# file = os.path.join("classifyre.json")
chemids = []

with open(file, "r") as f:
    ligands = json.load(f)
    pprint(len(ligands))
    # for ent in ligands:
    #     try:
    #         print(ent["identifier"])
    #         chemids.append(ent["identifier"])
    #     except:
    #         pprint(ent)

# pprint(chemids)