import json
import operator
from pprint import pprint
from time import time
from typing import Optional
import re
from typing import List
import warnings
from Bio import (
    BiopythonDeprecationWarning,
)
from fuzzysearch import find_near_matches
import numpy as np

from ribctl.lib.schema.types_ribosome import Polymer

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


#! Transposition methods

class SeqMatch:
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
        pprint(sourceseq)
        pprint(targetseq)

        _            = pairwise2.align.globalxx(self.src, self.tgt, one_alignment_only=True)
        self.src_aln = _[0].seqA
        self.tgt_aln = _[0].seqB
        #! The only thing that can happen hence is the insertion of gaps in the source sequence.

        self.aligned_ids = []

        for src_resid in self.src_ids:
            self.aligned_ids.append(self.forwards_match(self.src_aln, src_resid))

        self.aligned_ids = list(filter(lambda x: x != None, self.aligned_ids))

        for aln_resid in self.aligned_ids:
            if self.tgt_aln[aln_resid] == "-":
                continue
            self.tgt_ids.append(self.backwards_match(self.tgt_aln, aln_resid))

    def backwards_match(self, alntgt: str, resid: int):
        """Returns the target-sequence index of a residue in the (aligned) target sequence
        Basically, "count back ignoring gaps"
        """
        if resid > len(alntgt):
            raise IndexError(
                f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(alntgt)}"
            )
        counter_proper = 0
        for i, char in enumerate(alntgt):
            if i == resid:
                return counter_proper
            if char == "-":
                continue
            else:
                counter_proper += 1

    def forwards_match(self, alnsrc: str, resid: int):
        """Returns the index of a source-sequence residue in the aligned source sequence.
        Basically, "count forward including gaps"
        """

        count_proper = 0
        for alignment_indx, char in enumerate(alnsrc):
            if count_proper == resid:
                return alignment_indx
            if char == "-":
                continue
            else:
                count_proper += 1

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
    auth_seq_ids = {}

    for residue in nbr_residues:

        parent_chain = residue.get_parent()
        auth_asym_id = parent_chain.get_id()

        if auth_asym_id not in auth_seq_ids:
            auth_seq_ids[auth_asym_id] = [residue]
        else:
            auth_seq_ids[auth_asym_id].append(residue)

    RO = RibosomeOps(struct.get_id().upper())

    for chain_aaid, bound_residues in auth_seq_ids.items():
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


def create_residue_mapping_mask(canonical: str, structure: str) -> dict[int, int]:
    # Perform a global alignment
    alignments = pairwise2.align.globalxx(Seq(canonical), Seq(structure))
    aligned_canonical, aligned_structure = alignments[0][0], alignments[0][1]
    
    print('------------------------------------_***')
    pprint(aligned_canonical)
    pprint(aligned_structure)
    print('------------------------------------_***')
    # Convert aligned sequences to NumPy arrays
    can_array    = np.array(list(aligned_canonical))
    struct_array = np.array(list(aligned_structure))
    
    # Create indices for non-gap positions
    can_indices    = np.arange(len(can_array))[can_array != '-']
    struct_indices = np.cumsum(struct_array != '-') - 1
    
    # Create a mask for positions where both sequences have residues
    mask = (can_array != '-') 
    
    # Initialize the mapping dictionary with all canonical indices mapped to -1
    mapping = {i: -1 for i in range(len(canonical))}
    
    # Update the mapping for positions where both sequences have residues
    mapping.update(zip(can_indices[mask], struct_indices[mask]))
    
    return mapping




class SeqMap:

    mapping       : dict[int, int]

    seq_canonical : str
    seq_structural: str

    seq_canonical_aligned : str
    seq_structural_aligned: str


    def __init__(self, canonical: str, structure: str):
        self.seq_canonical  = canonical
        self.seq_structural = structure
        alignments          = pairwise2.align.globalxx(Seq(canonical), Seq(structure))

        aligned_canonical, aligned_structure = alignments[0][0], alignments[0][1]
        
        self.seq_canonical_aligned  = aligned_canonical
        self.seq_structural_aligned = aligned_structure

        mapping         = {}

        print("inspecting")
        print(self.seq_canonical_aligned)
        print(self.seq_structural_aligned)
        print(*zip(aligned_canonical, aligned_structure))
        # exit()
        
        canonical_index = 0
        structure_index = 0
        for canonical_char, structural_char in zip(aligned_canonical, aligned_structure):
            if canonical_char != '-':
                if structural_char != '-':
                    mapping[canonical_index] = structure_index
                    structure_index +=1  
                else:
                    mapping[canonical_index] = -1
                canonical_index +=1
            elif canonical_char == '-':
                warnings.warn(f"Unexpected gap in canonical sequence at aligned position {canonical_index}. This shouldn't happen with the original canonical sequence.")

        self.mapping = mapping

    def retrieve_index(self, key:int) ->int | None:
        "Get the STRUCTURAL sequence index corresponding to the CANONICAL sequence index <key> if any, otherwise None"
        if key not in self.mapping:
            raise KeyError(f"Key {key} not found in mapping")
        if self.mapping[key] == -1:
            return None
        return self.mapping[key]
        
    def retrieve_motif(self, keys:list[int]) -> tuple[str,str]:
        can_subseq    = ""
        struct_subseq = ""

        for i in keys:
            can_subseq = can_subseq + self.seq_canonical[i]
            struct_index = self.retrieve_index(i)
            if struct_index == None:
                struct_subseq = struct_subseq + "-"
            elif struct_index  != None:
                struct_subseq = struct_subseq + self.seq_structural[struct_index]
        if struct_subseq == "" or list(set(list(struct_subseq)))[0] == "-":
            raise ValueError("No structural sequence found for the given canonical sequence")
        return can_subseq, struct_subseq

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

#! Deperecated
def __bsite_transpose_motifs(
    source_rcsb_id: str,
    target_rcsb_id: str,
    binding_site: BindingSite,
    save: bool = False,
) -> LigandTransposition:

    start = time()

    def BiopythonChain_to_sequence(chain: Chain) -> tuple[str, dict[int, int]]:
        res: list[Residue] = [*chain.get_residues()]
        idx_auth_seq_id_map = {}
        seq = ""
        for idx, residue in enumerate(res):
            if residue.resname in [*AMINO_ACIDS.keys()]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
                idx_auth_seq_id_map[idx] = residue
            elif residue.resname in [*NUCLEOTIDES]:
                seq = seq + residue.resname
                idx_auth_seq_id_map[idx] = residue
            else:
                seq = seq + "-"
                idx_auth_seq_id_map[idx] = residue

        return seq, idx_auth_seq_id_map

    def extract_contiguous_motifs(
        bound_residues: list[tuple[int, str]]
    ) -> list[list[tuple[int, str]]]:
        bound_residues = sorted(
            bound_residues,
            key=lambda x: x[0],
        )
        motifs = []
        current_motif = []
        for i, (num, amino) in [*enumerate(bound_residues)]:
            if not current_motif:
                current_motif.append((num, amino))
            else:
                last_num = current_motif[-1][0]
                if num - last_num <= 3:  # Allow for up to 2 skipped numbers
                    current_motif.append((num, amino))
                else:
                    motifs.append(current_motif)
                    current_motif = [(num, amino)]

        if current_motif:
            motifs.append(current_motif)
        return motifs

    source_rcsb_id, target_rcsb_id = source_rcsb_id.upper(), target_rcsb_id.upper()
    source_polymers_by_poly_class = {}

    target_struct = RibosomeOps(target_rcsb_id).biopython_structure()
    source_struct = RibosomeOps(source_rcsb_id).biopython_structure()

    #! Source polymers
    for nbr_polymer in binding_site.chains:
        nbr_polymer = BindingSiteChain.model_validate(nbr_polymer)
        #! Skip if no nomenclature present ( can't do anything with it )
        if len(nbr_polymer.nomenclature) < 1:
            continue
        else:
            bound_residues_ids = [
                (resid, resname)
                for (resid, resname) in [
                    *map(
                        lambda x: (x.auth_seq_id, x.label_comp_id),
                        nbr_polymer.bound_residues,
                    )
                ]
            ]
            source_polymers_by_poly_class[nbr_polymer.nomenclature[0].value] = {
                # "seq"           : nbr_polymer.entity_poly_seq_one_letter_code_can,
                # "auth_asym_id"  : nbr_polymer.auth_asym_id,
                # "bound_residues": nbr_polymer.bound_residues,
                "polymer": nbr_polymer,
                #! Collect contiguous motifs for each polymer (to possibly seek them in the target)
                "motifs": extract_contiguous_motifs(bound_residues_ids),
            }

    #! Target polymers
    target_polymers: list[PredictedResiduesPolymer] = []
    for (
        nomenclature_class,
        source_polymer_with_motifs,
    ) in sorted(source_polymers_by_poly_class.items(), key=lambda x: x[0]):

        print("Processing chain [{}]".format(nomenclature_class))
        target_polymer: Polymer | None = RibosomeOps(
            target_rcsb_id
        ).get_poly_by_polyclass(nomenclature_class, 0)
        source_polymer: BindingSiteChain = source_polymers_by_poly_class[ nomenclature_class ]["polymer"]
        #! If no polymer of corresponding class is found, move on.
        if target_polymer == None:
            continue

        seq_src, idx_auth_map_src = BiopythonChain_to_sequence(
            source_struct[0][source_polymer.auth_asym_id]
        )
        seq_tgt, idx_auth_map_tgt = BiopythonChain_to_sequence(
            target_struct[0][target_polymer.auth_asym_id]
        )

        target_polymer_all_residues = []

        for i, motif in enumerate(source_polymer_with_motifs["motifs"]):

            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------
            motif_str = ""

            for _, residue_label in motif:
                motif_str = motif_str + ResidueSummary.three_letter_code_to_one( residue_label )

            if len(motif_str) <= 5:
                continue
            print("\tSource-motif {}: {}".format(i, motif_str))
            matches = find_near_matches(
                motif_str,
                seq_tgt,
                max_substitutions = 0,
                max_l_dist        = 1,
                max_insertions    = 2,
                max_deletions     = 0,
            )
            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------

            for j, match in enumerate(matches):
                target_motif_residues = []

                for i in range(match.start, match.end):
                    target_motif_residues.append(idx_auth_map_tgt[i])
                print( "\t\tTarget-motif match {}: {}".format( j, "".join( list( map( lambda x: ResidueSummary.three_letter_code_to_one( x.resname ), target_motif_residues)))))

                target_polymer_all_residues = [ *target_polymer_all_residues, *target_motif_residues, ]

        target_polymers.append(
            PredictedResiduesPolymer(
                polymer_class=nomenclature_class,
                source=PredictionSource(
                    source_seq=seq_src,
                    auth_asym_id=source_polymer.auth_asym_id,
                    source_bound_residues=[
                        ResidueSummary(
                            auth_seq_id=residue.auth_seq_id,
                            label_comp_id=residue.label_comp_id,
                            auth_asym_id=source_polymer.auth_asym_id,
                            label_seq_id=None,
                            full_id=None,
                            rcsb_id=source_rcsb_id,
                        )
                        for residue in source_polymer.bound_residues
                    ],
                ),
                target=PredictionTarget(
                    target_seq=seq_tgt,
                    auth_asym_id=target_polymer.auth_asym_id,
                    target_bound_residues=[
                        ResidueSummary(
                            auth_seq_id=residue.get_id()[1],
                            label_comp_id=residue.resname,
                            label_seq_id=None,
                            auth_asym_id=target_polymer.auth_asym_id,
                            full_id=None,
                            rcsb_id=target_rcsb_id,
                        )
                        for residue in target_polymer_all_residues
                    ],
                ),
            )
        )

    # ! at this point we have collected all the source polymers and corresponding target polymers.
    chains: list[BindingSiteChain] = []

    for c in target_polymers:
        poly = RibosomeOps(target_rcsb_id).get_poly_by_auth_asym_id(
            c.target.auth_asym_id
        )
        if poly == None:
            raise ValueError("Polymer not found in target structure")
        bsite_chain = BindingSiteChain(
            **poly.model_dump(), bound_residues=c.target.target_bound_residues
        )
        chains.append(bsite_chain)

    end = time()

    _ = LigandTransposition(
        constituent_chains=target_polymers,  # this is the info about each individual pair of polymers manipulations
        source=source_rcsb_id,
        target=target_rcsb_id,
        purported_binding_site=BindingSite(  # this is the result, the datastructure that gets sent the fronted
            chains=chains,
            ligand=binding_site.ligand,
            radius=binding_site.radius,
            source=binding_site.source,
        ),
    )
    # pprint(_)
    print("Elapsed time: ", end - start)
    return _

def bsite_transpose(
  source_rcsb_id        : str,
  target_rcsb_id        : str,
  binding_site          : BindingSite,
  save                  : bool = False,
) -> LigandTransposition: 

    start = time()

    def BiopythonChain_to_sequence(chain: Chain) -> str:
        res: list[Residue] = [*chain.get_residues()]
        seq = ""
        for idx, residue in enumerate(res):
            if residue.resname in [*AMINO_ACIDS.keys()]:
                seq = seq + ResidueSummary.three_letter_code_to_one(residue.resname)
            elif residue.resname in [*NUCLEOTIDES]:
                seq = seq + residue.resname
            else:
                continue
                seq = seq + "-"

        return seq

    source_rcsb_id, target_rcsb_id = source_rcsb_id.upper(), target_rcsb_id.upper()

    source_polymers_by_poly_class = {}

    target_struct = RibosomeOps(target_rcsb_id).biopython_structure()
    source_struct = RibosomeOps(source_rcsb_id).biopython_structure()

    target_profile = RibosomeOps(target_rcsb_id).profile()
    source_profile = RibosomeOps(source_rcsb_id).profile()

    #* Work out a mapping between the structural and the canonical sequences






    #* Work out a mapping between the structural and the canonical sequences
    #! Source polymers
    for nbr_polymer in binding_site.chains:
        if "uS12" not in nbr_polymer.nomenclature :
            continue
        nbr_polymer = BindingSiteChain.model_validate(nbr_polymer)
        #! Skip if no nomenclature present ( can't do anything with it )
        if len(nbr_polymer.nomenclature) < 1:
            continue
        else:
            # ! Bound residues  [in STRUCTURE SPACE]
            # bound_residues_ids = [ (resid, resname) for (resid, resname) in [ *map( lambda x: (x.auth_seq_id, x.label_comp_id), nbr_polymer.bound_residues)]]
            canonical_seq_src  = RibosomeOps(source_rcsb_id).get_poly_by_auth_asym_id(nbr_polymer.auth_asym_id).entity_poly_seq_one_letter_code_can
            structural_seq_src = BiopythonChain_to_sequence(source_struct[0][nbr_polymer.auth_asym_id] )

            # pprint(structural_seq_src[:40])
            M = SeqMap( canonical_seq_src, structural_seq_src)
            print("\n\nCanonical sequence")
            pprint(M.seq_canonical)
            print("Structural sequence")
            pprint(M.seq_structural)
            print("\n Both, aligned:")
            print(M.seq_canonical_aligned)
            print(M.seq_structural_aligned)
            print("\t\t\t\t*************")

            sub_ixs       = [120, 50, 22, 33, 44]
            can_subseq    = ""
            struct_subseq = ""

            pprint(M.mapping)

            can, stru = M.retrieve_motif(sub_ixs)
            print("Canonical subsequence", )
            print(can)
            print("structu")
            print(stru)


            source_polymers_by_poly_class[nbr_polymer.nomenclature[0].value] = {
                # "seq"           : nbr_polymer.entity_poly_seq_one_letter_code_can,
                # "auth_asym_id"  : nbr_polymer.auth_asym_id,
                # "bound_residues": nbr_polymer.bound_residues,
                "polymer": nbr_polymer,
                #! Collect contiguous motifs for each polymer (to possibly seek them in the target)
                # "motifs": extract_contiguous_motifs(bound_residues_ids),
            }

    exit()
    #! Target polymers
    target_polymers: list[PredictedResiduesPolymer] = []
    for ( nomenclature_class, source_polymer_with_motifs, ) in sorted(source_polymers_by_poly_class.items(), key=lambda x: x[0]):
        print("Processing chain [{}]".format(nomenclature_class))
        target_polymer: Polymer | None = RibosomeOps( target_rcsb_id ).get_poly_by_polyclass(nomenclature_class, 0)
        source_polymer: BindingSiteChain = source_polymers_by_poly_class[ nomenclature_class ]["polymer"]
        #! If no polymer of corresponding class is found, move on.
        if target_polymer == None:
            continue

        structural_seq_src, idx_auth_map_src = BiopythonChain_to_sequence( source_struct[0][source_polymer.auth_asym_id] )
        seq_tgt, idx_auth_map_tgt = BiopythonChain_to_sequence( target_struct[0][target_polymer.auth_asym_id] )

        target_polymer_all_residues = []

        for i, motif in enumerate(source_polymer_with_motifs["motifs"]):

            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------
            motif_str = ""
            for _, residue_label in motif:
                motif_str = motif_str + ResidueSummary.three_letter_code_to_one( residue_label )

            if len(motif_str) <= 5:
                continue
            print("\tSource-motif {}: {}".format(i, motif_str))
            matches = find_near_matches(
                motif_str,
                seq_tgt,
                max_substitutions = 0,
                max_l_dist        = 1,
                max_insertions    = 2,
                max_deletions     = 0,
            )
            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------

            for j, match in enumerate(matches):
                target_motif_residues = []

                for i in range(match.start, match.end):
                    target_motif_residues.append(idx_auth_map_tgt[i])
                print( "\t\tTarget-motif match {}: {}".format( j, "".join( list( map( lambda x: ResidueSummary.three_letter_code_to_one( x.resname ), target_motif_residues)))))

                target_polymer_all_residues = [ *target_polymer_all_residues, *target_motif_residues, ]

        target_polymers.append(
            PredictedResiduesPolymer(
                polymer_class=nomenclature_class,
                source=PredictionSource(
                    source_seq=structural_seq_src,
                    auth_asym_id=source_polymer.auth_asym_id,
                    source_bound_residues=[
                        ResidueSummary(
                            auth_seq_id=residue.auth_seq_id,
                            label_comp_id=residue.label_comp_id,
                            auth_asym_id=source_polymer.auth_asym_id,
                            label_seq_id=None,
                            full_id=None,
                            rcsb_id=source_rcsb_id,
                        )
                        for residue in source_polymer.bound_residues
                    ],
                ),
                target=PredictionTarget(
                    target_seq=seq_tgt,
                    auth_asym_id=target_polymer.auth_asym_id,
                    target_bound_residues=[
                        ResidueSummary(
                            auth_seq_id=residue.get_id()[1],
                            label_comp_id=residue.resname,
                            label_seq_id=None,
                            auth_asym_id=target_polymer.auth_asym_id,
                            full_id=None,
                            rcsb_id=target_rcsb_id,
                        )
                        for residue in target_polymer_all_residues
                    ],
                ),
            )
        )

    # ! at this point we have collected all the source polymers and corresponding target polymers.
    chains: list[BindingSiteChain] = []

    for c in target_polymers:
        poly = RibosomeOps(target_rcsb_id).get_poly_by_auth_asym_id( c.target.auth_asym_id )
        if poly == None:
            raise ValueError("Polymer not found in target structure")
        bsite_chain = BindingSiteChain( **poly.model_dump(), bound_residues=c.target.target_bound_residues )
        chains.append(bsite_chain)

    end = time()

    _ = LigandTransposition(
        constituent_chains=target_polymers,  # this is the info about each individual pair of polymers manipulations
        source=source_rcsb_id,
        target=target_rcsb_id,
        purported_binding_site=BindingSite(  # this is the result, the datastructure that gets sent the fronted
            chains=chains,
            ligand=binding_site.ligand,
            radius=binding_site.radius,
            source=binding_site.source,
        ),
    )
    # pprint(_)
    print("Elapsed time: ", end - start)
    return _


# Old alignment/matching logic:

# for nomenclature_class, _ in source_polymers_by_poly_class.items():
#     if nomenclature_class not in target_polymers_by_poly_class:
#         continue

# src_ids = source_polymers_by_poly_class[nomenclature_class]["ids"]
# src     = source_polymers_by_poly_class[nomenclature_class]["seq"]

# tgt     = target_polymers_by_poly_class[nomenclature_class]["seq"]

# src_auth_asym_id = source_polymers_by_poly_class[nomenclature_class][ "auth_asym_id" ]
# tgt_auth_asym_id = target_polymers_by_poly_class[nomenclature_class][ "auth_asym_id" ]

# sq               = SeqMatch(src, tgt, src_ids)

# src_aln = ( sq.src_aln )  # <--- aligned source      sequence (with                        gaps)
# tgt_aln = ( sq.tgt_aln )  # <--- aligned tgt         sequence (with                        gaps)
# aln_ids = ( sq.aligned_ids )  # <--- ids     corrected   for                                   gaps
# tgt_ids = ( sq.tgt_ids )  # <--- ids     backtracted to the target polymer (accounting for gaps)
