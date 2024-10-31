from pprint import pprint
import sys
from typing import NewType, TypeVar
from Bio.PDB.Residue import Residue
from ribctl.etl.etl_assets_ops import RibosomeOps
from ribctl.lib.libbsite import map_motifs, bsite_ligand, bsite_transpose
from ribctl.lib.libseq import SequenceMappingContainer
from ribctl.lib.schema.types_binding_site import (
    PredictionTarget,
    ResiduesMapping,
    ResidueSummary,
)
from ribctl.lib.schema.types_ribosome import PolymerClass
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.lib.libmsa import Fasta
from concurrent.futures import ProcessPoolExecutor
from multiprocessing import Pool
from functools import partial

# ? Ligand Projections
# Ok let's not complicate things. The prototype is thus:
# take a motif in 7K00 protein, say uL10

# 1. To Record:
# - locus name and its anchor structure/chain/residue
# - project anchor into the MSA of homologs of other structure
# - project anchor into each other structure PAIRWISE

# 2. Backtrack to structural files indices for each structure, record as a row in a table.

def process_single_structure(
    structure_pair: tuple[str,str],
    target_rcsb_id: str  = "7K00",
    radius: float = 15
) -> list[tuple[PolymerClass, PredictionTarget]]:
    """Process a single structure and return list of (polymer_class, target) pairs"""

    source_rcsb_id, chem_id = structure_pair
    print(f"Mapping {source_rcsb_id}/{chem_id}(source) into {target_rcsb_id} (target)")

    bsite_source = bsite_ligand(chem_id, source_rcsb_id, radius)
    bsite_target = bsite_transpose(source_rcsb_id, target_rcsb_id, bsite_source)

    results = []
    for chain in bsite_target.constituent_chains:
        results.append((chain.polymer_class, chain.target))
    return results

def prepare_mapping_sources(
    source_structures: list[tuple[str, str]], n_processes: int = 8
) -> dict[PolymerClass, list[PredictionTarget]]:
    per_class_registry: dict[PolymerClass, list[PredictionTarget]] = {}
    with Pool(processes = n_processes) as pool:
        all_results = pool.map(process_single_structure, source_structures)

    for structure_results in all_results:
        for polymer_class, target in structure_results:
            if polymer_class not in per_class_registry:
                per_class_registry[polymer_class] = [target]
            else:
                per_class_registry[polymer_class].append(target)
    return per_class_registry

def compact_class(
    projections: list[PredictionTarget],
    number_of_sources: int,
) -> dict[int, float]:

    weights = {}

    for mapping in projections:
        residues = mapping.target_bound_residues
        for residue in residues:
            if residue.auth_seq_id not in weights:
                weights.update({residue.auth_seq_id: 1})
            else:
                weights[residue.auth_seq_id] += 1

    [ weights.update({x:y})  for x,y in list(map(lambda item: (item[0], item[1] / number_of_sources), weights.items())) ] 
    return weights




# class MotifsMapManyToOne:
#     """Given multiple source sequences and residue ranges within them, project the ranges onto a single target sequence"""

#     target_chain : BiopythonChain
#     target_rcsb_id: RCSB_ID
#     polymer_class: PolymerClass
#     projections  : list[PredictionTarget]


#     def __init__(
#         self,
#       target_chain  : BiopythonChain,
#       target_rcsb_id: RCSB_ID,
#       motif_mappings: list[PredictionTarget],
#       polymer_class : PolymerClass
#     ) -> None       :
#         self.target_chain   = target_chain
#         self.target_rcsb_id = target_rcsb_id
#         self.polymer_class  = PolymerClass(polymer_class)
#         self.projections    = motif_mappings

#     def get_weights(self, number_of_instances) -> dict[int, list[RCSB_ID]]:

#         return {}

# return MotifsMapManyToOne(
#     target_chain,
#     target_rcsb_id,
#     polymer_class,


# )
# target chain

# ...

# target_chain_aaid = RibosomeOps(target_rcsb_id).get_poly_by_polyclass(PolymerClass("16SrRNA")).auth_asym_id
# target_chain = RibosomeOps(target_rcsb_id).biopython_structure()[0][target_chain_aaid]