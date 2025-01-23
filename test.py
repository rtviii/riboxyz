# Create Fasta object with sequences
import json
import os
from pprint import pprint
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from neo4j_ribosome.db_lib_reader import Neo4jReader
from ribctl import ASSETS, ASSETS_PATH
from ribctl.lib.libbsite import bsite_transpose, get_lig_bsite
from ribctl.lib.libmsa import Fasta, FastaBalancer
from ribctl.lib.libseq import (
    SequenceMappingContainer,
    get_conservation_scores,
    map_alignment_to_structure,
)
from ribctl.lib.libtax import Taxid, get_lineage_distance, print_lineage_comparison
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


# file = os.path.join(ASSETS_PATH, "ligands", "ligand_classes.json")
# # file = os.path.join("classifyre.json")
# chemids = []

# with open(file, "r") as f:
#     ligands = json.load(f)
#     pprint(len(ligands))
#     # for ent in ligands:
#     #     try:
#     #         print(ent["identifier"])
#     #         chemids.append(ent["identifier"])
#     #     except:
#     #         pprint(ent)

# # pprint(chemids)



# bsite = get_lig_bsite('YQM', RibosomeOps('7M4W').assets.biopython_structure(), 10.0)
# res = bsite_transpose('7M4W','7K00',bsite)
all_ligand_sources = Neo4jReader(Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)).list_ligands()
# pprint(all_ligand_ids)
RO = RibosomeOps('7K00')

for lig_source in all_ligand_sources:
    distance         = 999999
    closest_to_ecoli = None
    for elem in lig_source[1:]:
        for ss in elem:
            distance = get_lineage_distance(RO.taxid,ss['tax_node']['ncbi_tax_id'])

    # exit()

    # pprint(s)
    # Taxid.relative_distances_to_taxid(RO.
    
    # s['ncbi_tax_id']
    # exit()
# for CHEMID in all_ligand_ids:
#     bsite = get_lig_bsite(CHEMID, RibosomeOps('7M4W').assets.biopython_structure(), 10.0)
#     res   = bsite_transpose('7M4W','7K00',bsite)