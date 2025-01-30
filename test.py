# Create Fasta object with sequences
import json
import os
from pprint import pprint
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from neo4j_ribosome.db_lib_reader import Neo4jReader
from ribctl import ASSETS, ASSETS_PATH
from ribctl.global_ops import GlobalOps
from ribctl.lib.libbsite import bsite_transpose, extract_ligand_to_mmcif, get_lig_bsite
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


# lig_chemid = 'SAH'
# source_id  = '6SG9'
# tartget_id = '7K00'

# bsite = get_lig_bsite(lig_chemid, RibosomeOps(source_id).assets.biopython_structure(), 10.0)
# pprint(bsite)
# res   = bsite_transpose(source_id,tartget_id,bsite)
# print(res)

# all_ligand_sources = Neo4jReader(Neo4jAdapter(NEO4J_URI, NEO4J_USER, NEO4J_CURRENTDB)).list_ligands()
# RO                 = RibosomeOps('7K00')

# for lig_source in all_ligand_sources:
#     distance_0         = 999999
#     closest_to_ecoli = None

#     for elem in lig_source[1:]:
#         for ss in elem:
#             distance = get_lineage_distance(RO.taxid,ss['tax_node']['ncbi_tax_id'])
#             if distance <  distance_0:
#                 distance_0       = distance
#                 closest_to_ecoli = ss

#     # pprint(closest_to_ecoli)
#     chemical_id = lig_source[0]['chemicalId']
#     source_id   = closest_to_ecoli['rcsb_id']
#     outpath     = RO.assets.paths.binding_site_prediction(chemical_id, source_id)

#     if outpath is not None and os.path.exists(outpath):
#         print(f"File exists: {outpath}")
#     else:
#         print(f"Closest to Ecoli for ligand {lig_source[0]['chemicalId']}: {closest_to_ecoli['rcsb_id']} {closest_to_ecoli['tax_node']['ncbi_tax_id']} {closest_to_ecoli['tax_node']['scientific_name']} {distance_0}")
#         bsite = get_lig_bsite(chemical_id, RibosomeOps(source_id).assets.biopython_structure(), 10.0)
#         res   = bsite_transpose(source_id,'7K00',bsite)
#         with open(outpath, "w") as f:
#             json.dump(res.model_dump(), f)
#             print("Saved: ", outpath)

from neo4j_ribosome.db_lib_reader import dbqueries

print(dbqueries.ligands_per_structure())