import json
import operator
from time import time
from typing import Optional
import re
from typing import List
import warnings
from Bio import (
    BiopythonDeprecationWarning,
)
from fuzzysearch import find_near_matches
import pandas as pd
from ribctl.lib.schema.types_ribosome import Polymer
sys.path.append("/home/rtviii/dev/riboxyz")
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import pairwise2
from ribctl.lib.schema.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
    BindingSiteChain,
    LigandTransposition,
    PredictedResiduesPolymer,
    PredictionAlignments,
    PredictionSource,
    PredictionTarget,
)
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.Seq import Seq
from Bio.PDB.Structure import Structure
from ribctl.etl.etl_assets_ops import RibosomeOps
from ribctl.lib.schema.types_binding_site import (
    BindingSite,
    BindingSiteChain,
    ResidueSummary,
)
import os
import sys
from collections import defaultdict
from ribctl import ASSETS_PATH


def lig_get_chemical_categories():
    ligands_classification_path = os.path.join(ASSETS_PATH, "ligands", "ligand_chemical_categories.csv")
    df           = pd.read_csv(ligands_classification_path)
    reverse_dict = defaultdict(list)
    for _, row in df.iterrows():
        ligand = row['Ligand']
        category = row['Category']
        reverse_dict[category].append(ligand)
    
    return dict(reverse_dict)


#! Transposition methods
class BiopythonChain(Chain):
    chain:Chain 
    flat_index_to_residue_map    : dict[int, Residue]
    auth_seq_id_to_flat_index_map: dict[int, int]

    def __init__(self, chain:Chain):
        self.chain = chain

    @property
    def primary_sequence(self) -> tuple[str, dict]:

        represent_noncanonical_as:Optional[str]="."
        seq = ""
        auth_seq_id_to_primary_ix = {}
        for ix,residue in enumerate( self.chain.get_residues() ):
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

        flat_index_to_residue_map     = {}
        auth_seq_id_to_flat_index_map = {}
        seq                           = ""
        flat_index                         = 0
        for residue in res:
            if residue.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]:
                seq                                                = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                flat_index_to_residue_map[flat_index]              = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = flat_index
                flat_index +=1
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
        self.tgt    : str       = targetseq
        self.tgt_ids: list[int] = []

        _ = pairwise2.align.globalxx(self.src, self.tgt, one_alignment_only=True)

        self.src_aln = _[0].seqA
        self.tgt_aln = _[0].seqB


        self.aligned_ids = []

        for src_resid in self.src_ids:
            self.aligned_ids.append(self.forwards_match(self.src_aln, src_resid) )

        for aln_resid in self.aligned_ids:
            tgt_aln_index = self.backwards_match(self.tgt_aln, aln_resid)
            if tgt_aln_index == None:
                continue
            else:
                self.tgt_ids.append(tgt_aln_index)
        
        print("[Source Aligned]\t",self.hl_ixs(self.src_aln, ixs=self.aligned_ids))
        print("[Target Aligned]\t",self.hl_ixs(self.tgt_aln, ixs=self.aligned_ids))

        

    def forwards_match(self, aligned_source_sequence: str, original_residue_index: int)->int:
        """Returns the index of a source-sequence residue in the aligned source sequence. Basically, "count forward including gaps"
        """
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

        raise ValueError(f"Residue with index {original_residue_index} not found in the aligned source sequence after full search. Logical errory, likely.")

    def backwards_match(self, aligned_target_sequence: str, aligned_residue_index: int)->int|None:
        """Returns the target-sequence index of a residue in the [aligned] target sequence. Basically, "count back ignoring gaps"
        """
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
        self.seq_canonical = canonical
        self.seq_structural = structure
        alignments = pairwise2.align.globalxx(Seq(canonical), Seq(structure))

        aligned_canonical, aligned_structure = alignments[0][0], alignments[0][1]

        self.seq_canonical_aligned = aligned_canonical
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
        can_subseq = ""
        struct_subseq = ""

        for i in keys:
            can_subseq = can_subseq + self.seq_canonical[i]
            struct_index = self.retrieve_index(i)
            if struct_index == None:
                struct_subseq = struct_subseq + "-"
            elif struct_index != None:
                struct_subseq = struct_subseq + self.seq_structural[struct_index]
        if struct_subseq == "" or list(set(list(struct_subseq)))[0] == "-":
            raise ValueError(
                "No structural sequence found for the given canonical sequence"
            )
        return can_subseq, struct_subseq

def get_lig_bsite(
    lig_chemid: str,
    struct: Structure,
    radius: float,
) -> BindingSite:
    """KDTree search the neighbors of a given list of residues (which constitue a ligand)
    and return unique
    """
    # Make sure only the first assembly is used if multiple are in the file.
    md: Model = [*struct.get_models()][0]
    assemblies = (
        RibosomeOps(struct.get_id().upper()).profile().get_polymers_by_assembly()
    )
    # If there are two or more assemblies, delete chains belonging to all but the first one.
    if len(assemblies.items()) > 1:
        for i in range(len(assemblies.items()) - 1):
            for chain_aaid in [*assemblies.items()][i + 1][1]:
                md.detach_child(chain_aaid)

    ns = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []
    ligand_residues = list(
        filter(lambda x: x.get_resname() == lig_chemid, list(struct.get_residues()))
    )

    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
            found_nbrs = ns.search(atom.get_coord(), radius, level="R")
            nbr_residues.extend(found_nbrs)

    nbr_residues: list[Residue] = list(set(nbr_residues))
    nbr_chains = []

    # TODO: grab the PROXIMATE idx by `enumerate`ing the chain residues
    # It doesn't fucking matter as long as you are unable to match residue comp ID 1-to-1
    # TODO: Just try constructing the sequence from the biopython one and aligning (with pairwise2) into the canonical one.
    # TODO:          (then just subtract the gaps from the canonical one)
    # auth_seq_id_to_label_seq_id_mapping_by_chain = {}
    nbr_residues_by_chain_aaid = {}

    for residue in nbr_residues:

        parent_chain = residue.get_parent()
        auth_asym_id = parent_chain.get_id()
        if auth_asym_id not in nbr_residues_by_chain_aaid:
            nbr_residues_by_chain_aaid[auth_asym_id] = [residue]
        else:
            nbr_residues_by_chain_aaid[auth_asym_id].append(residue)

    RO = RibosomeOps(struct.get_id().upper())

    # pprint(sorted( nbr_residues_by_chain_aaid['L'], key=lambda x: x.get_id()[1] ))
    for chain_aaid, bound_residues in nbr_residues_by_chain_aaid.items():
        polymer = RO.get_poly_by_auth_asym_id(chain_aaid)
        if polymer == None:
            raise ValueError(
                f"Polymer with auth_asym_id {chain_aaid} not found in structure. Logic error."
            )
        bound_residues = sorted(
            [
                ResidueSummary(
                    full_id=None,
                    auth_asym_id=chain_aaid,
                    label_comp_id=residue.resname,
                    auth_seq_id=residue.get_id()[1],
                    label_seq_id=None,
                    rcsb_id=struct.get_id().upper(),
                )
                for residue in bound_residues
            ],
            key=operator.attrgetter("auth_seq_id"),
        )
        nbr_chains.append(
            BindingSiteChain(**polymer.model_dump(), bound_residues=bound_residues)
        )

    return BindingSite(
        chains=nbr_chains,
        ligand=lig_chemid,
        radius=radius,
        source=struct.get_id().upper(),
    )

def bsite_ligand(
    chemicalId: str, rcsb_id: str, radius: float, save: bool = False
) -> BindingSite:

    chemicalId = chemicalId.upper()
    _structure_cif_handle = RibosomeOps(rcsb_id).biopython_structure()
    binding_site_ligand = get_lig_bsite(chemicalId, _structure_cif_handle, radius)

    if save:
        with open(RibosomeOps(rcsb_id).paths.binding_site(chemicalId), "w") as f:
            json.dump(binding_site_ligand.model_dump(), f)

    return binding_site_ligand

def bsite_transpose(
    source_rcsb_id: str,
    target_rcsb_id: str,
    binding_site: BindingSite,
    save: bool = False,
) -> LigandTransposition:

    start = time()

    def BiopythonChain_to_sequence(chain: Chain) -> tuple[str, dict, dict]:
        res: list[Residue] = [*chain.get_residues()]
        
        flat_index_to_residue_map     = {}
        auth_seq_id_to_flat_index_map = {}
        seq                           = ""
        index                         = 0

        for residue in res:
            if residue.resname in [*AMINO_ACIDS.keys()]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                index += 1
                flat_index_to_residue_map[index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = index
            elif residue.resname in [*NUCLEOTIDES]:
                seq = seq + residue.resname
                index += 1
                flat_index_to_residue_map[index] = residue
                auth_seq_id_to_flat_index_map[residue.get_id()[1]] = index
            else:
                continue
                seq = seq + "-"

        return seq, flat_index_to_residue_map, auth_seq_id_to_flat_index_map

    source_rcsb_id, target_rcsb_id = source_rcsb_id.upper(), target_rcsb_id.upper()
    source_polymers_by_poly_class = {}

    target_struct = RibosomeOps(target_rcsb_id).biopython_structure()
    source_struct = RibosomeOps(source_rcsb_id).biopython_structure()

    source_ops = RibosomeOps(source_rcsb_id)
    source_profile = source_ops.profile()

    target_ops = RibosomeOps(target_rcsb_id)
    target_profile = target_ops.profile()

    # * Work out a mapping between the structural and the canonical sequences

    # * Work out a mapping between the structural and the canonical sequences
    #! Source polymers

    chain_mappings      = []
    bindign_site_chains = []

    for nbr_polymer in binding_site.chains:
        nbr_polymer = BindingSiteChain.model_validate(nbr_polymer)
        #! Skip if no nomenclature present ( can't do anything with it )
        if len(nbr_polymer.nomenclature) < 1:
            continue
        target_polymer = target_ops.get_chain_by_polymer_class( nbr_polymer.nomenclature[0], 0 )
        if target_polymer == None:
            continue
        
        print("\n\n\t\t [{}] ".format(nbr_polymer.nomenclature[0]))
        bpchain_source = BiopythonChain(source_struct[0][nbr_polymer.auth_asym_id])
        bpchain_target = BiopythonChain(target_struct[0][target_polymer.auth_asym_id])
        #! SOURCE MAPS
        [
            src_flat_structural_seq,
            src_flat_idx_to_residue_map,
            src_auth_seq_id_to_flat_index_map,
        ] = bpchain_source.flat_sequence
        #! TARGET MAPS
        [ 
            tgt_flat_structural_seq,
            tgt_flat_idx_to_residue_map,
            tgt_auth_seq_id_to_flat_index_map
        ] = bpchain_target.flat_sequence

        # ! Bound residues  [in STRUCTURE SPACE]
        src_bound_auth_seq_idx = [ (residue.auth_seq_id, residue.label_comp_id) for residue in nbr_polymer.bound_residues ]
        # ! Bound residues  [in STRUCTURE SPACE]
        
        primary_seq_source, auth_seq_to_primary_ix_source = bpchain_source.primary_sequence
        primary_seq_target, auth_seq_to_primary_ix_target = bpchain_target.primary_sequence
        print("[Source Primary]\t",SeqPairwise.hl_ixs(primary_seq_source, [ auth_seq_to_primary_ix_source[index] for index, label in src_bound_auth_seq_idx]))
        src_bound_flat_indices = [ src_auth_seq_id_to_flat_index_map[index] for index, label in filter( lambda x: x[1] in [*NUCLEOTIDES, *AMINO_ACIDS.keys()],src_bound_auth_seq_idx)]
        print("[Source Flat   ]\t",SeqPairwise.hl_ixs(src_flat_structural_seq, src_bound_flat_indices))
        print("- - - - Alignment- - - ")
        M                      = SeqPairwise(src_flat_structural_seq, tgt_flat_structural_seq, src_bound_flat_indices)
        print("- - - - Alignment- - - ")
        tgt_bound_flat_indices = M.tgt_ids
        tgt_bound_residues     = [ tgt_flat_idx_to_residue_map[idx] for idx in tgt_bound_flat_indices ]
        print("[Target Flat   ]\t", SeqPairwise.hl_ixs(tgt_flat_structural_seq, tgt_bound_flat_indices))
        print("[Target Primary]\t", SeqPairwise.hl_ixs(primary_seq_target,[auth_seq_to_primary_ix_target[residue.get_id()[1]] for residue in tgt_bound_residues] ))


        polymer_pair =PredictedResiduesPolymer(
            polymer_class=nbr_polymer.nomenclature[0],
            source=PredictionSource(
                    source_seq            = primary_seq_source,
                    auth_asym_id          = nbr_polymer.auth_asym_id,
                    source_bound_residues = [
                        ResidueSummary(
                            auth_seq_id   = residue.auth_seq_id,
                            label_comp_id = residue.label_comp_id,
                            auth_asym_id  = nbr_polymer.auth_asym_id,
                            label_seq_id  = None,
                            full_id       = None,
                            rcsb_id       = source_rcsb_id,
                        )
                        for residue in nbr_polymer.bound_residues
                    ],
                ),
                target=PredictionTarget(
                    target_seq            = primary_seq_target,
                    auth_asym_id          = target_polymer.auth_asym_id,
                    target_bound_residues = [
                        ResidueSummary(
                            auth_seq_id   = residue.get_id()[1],
                            label_comp_id = residue.resname,
                            label_seq_id  = None,
                            auth_asym_id  = target_polymer.auth_asym_id,
                            full_id       = None,
                            rcsb_id       = target_rcsb_id,
                        )
                        for residue in tgt_bound_residues
                    ],
                )

        )
        chain_mappings.append(polymer_pair)
        bindign_site_chains.append(BindingSiteChain(**target_polymer.model_dump(), bound_residues=[
                        ResidueSummary(
                            auth_seq_id   = residue.get_id()[1],
                            label_comp_id = residue.resname,
                            label_seq_id  = None,
                            auth_asym_id  = target_polymer.auth_asym_id,
                            full_id       = None,
                            rcsb_id       = target_rcsb_id,
                        )
                        for residue in tgt_bound_residues
                    ]))

    _ = LigandTransposition(
        source=source_rcsb_id,
        target=target_rcsb_id,
        constituent_chains=chain_mappings,  # this is the info about each individual pair of polymers manipulations
        purported_binding_site=BindingSite(  # this is the result, the datastructure that gets sent the fronted
            chains=bindign_site_chains,
            ligand=binding_site.ligand,
            radius=binding_site.radius,
            source=binding_site.source,
        ),
    )
    return _


# -----
# Projection workflows 
# The idea is to highlight the residues in the target structure based on a number of extant
#
# - we shouldn't reuse the same ligand twice (if there are multiple ligand instances available across structures -- select the closest taxonomically)
# Thinking about the ligands now. Last time we spoke, i think we landed on this idea of a "global" ligand chart i.e. given a random structure --  point out "hot pockets" given all the ligands we already have on hand. Did i get that right?
# How do i present that in practice?  Given that classes of ligands bind to ~the same locus ptc/tunnel lower/tunnel upper/decoding) i'll combine these predictions into "classes". The user would be able to select an individual ligand as they can now, but now also a class of ligands simultaneously from multiple structures, say "aminoglycosides" or "tetracyclines" and so on.
# It's not suuper taxing to classify the stuff we have, someone has already done the work : https://jcheminf.biomedcentral.com/articles/10.1186/s13321-016-0174-y
# Optics of it aside, i want to press you for why this is a potentially useful feature to have? I'll do it for sure, I think it is a very nice idea.
# I'm just wondering -- what's the next step to make it practical, open a door for people to build on it? Concrete example: we have 10 antibiotics in the class "Aminoglycosides" across 10 structures of 5-7 different species . I do the msa, mapping stuff. Now 10 extant (source) binding pockets are mapped into 1 target sequence. There is some overlap between source pockets, but also some spread in terms of what residues they get mapped into. How to "export" that? Does it tell anyone anything useful? Which ones are more relevant for drug design? (edited) 
# -----


def retrieve_ligands():

