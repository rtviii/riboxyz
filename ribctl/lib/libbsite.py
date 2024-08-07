import json
import operator
from pprint import pprint
from typing import Optional
import re
from typing import List
import warnings
from Bio import (
    BiopythonDeprecationWarning,
)
from fuzzysearch import find_near_matches

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
        """A container for origin and target sequences when matching the resiudes of a ligand binding site
        to another protein's sequence through BioSeq's Align
        """

        # *  indices of the ligand-facing residues in the source sequence.
        self.src: str = sourceseq
        self.src_ids: list[int] = source_residues

        # * Indices of predicted residues in target sequence. To be filled.
        self.tgt: str = targetseq
        self.tgt_ids: list[int] = []

        _ = pairwise2.align.globalxx(self.src, self.tgt, one_alignment_only=True)
        self.src_aln = _[0].seqA
        self.tgt_aln = _[0].seqB

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
                    parent_auth_asym_id=chain_aaid,
                    resname=residue.resname,
                    auth_seq_id=residue.get_id()[1],
                    label_seq_id=None,
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
        if len(nbr_polymer.nomenclature) < 1:
            continue
        else:
            biopython_chain: Chain = source_struct[0][nbr_polymer.auth_asym_id]
            bound_residues_ids = [
                (resid, resname)
                for (resid, resname) in [
                    *map(
                        lambda x: (x.auth_seq_id, x.resname), nbr_polymer.bound_residues
                    )
                ]
            ]
            source_polymers_by_poly_class[nbr_polymer.nomenclature[0].value] = {
                "seq": nbr_polymer.entity_poly_seq_one_letter_code_can,
                "auth_asym_id": nbr_polymer.auth_asym_id,
                "bound_residues": bound_residues_ids,
                "motifs": extract_contiguous_motifs(bound_residues_ids),
            }

    #! Target polymers
    target_polymers_by_poly_class = {}
    predicted_chains: list[PredictedResiduesPolymer] = []

    for nomenclature_class, nbr_polymer in source_polymers_by_poly_class.items():
        target_polymer = RibosomeOps(target_rcsb_id).get_poly_by_polyclass(
            nomenclature_class, 0
        )
        if target_polymer == None:
            continue
        # print(
        #     "\n\nMatched source chain {}.{} to target_polymer.auth_asym_id {}".format(
        #         nomenclature_class,
        #         source_polymers_by_poly_class[nomenclature_class]["auth_asym_id"],
        #         target_polymer.auth_asym_id,
        #     )
        # )

        seq_src, idx_auth_map_src = BiopythonChain_to_sequence(
            source_struct[0][nbr_polymer["auth_asym_id"]]
        )
        seq_tgt, idx_auth_map_tgt = BiopythonChain_to_sequence(
            target_struct[0][target_polymer.auth_asym_id]
        )

        all_motifs_residues = []
        for motif in nbr_polymer["motifs"]:
            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------
            motif_str = "".join( [ResidueSummary.three_letter_code_to_one(amino) for _, amino in motif] )
            if len(motif_str) <= 5:
                continue
            matches = find_near_matches( motif_str, seq_tgt, max_substitutions=0, max_l_dist=1, max_insertions=2, max_deletions=0, )
            # ! -------------------------------------------- SEARCH PARAMS --------------------------------------------------

            # TODO : Use all matches.
            # TODO : Don't discard matches. The clustering should take care of this.
            if len(matches) != 1:
                continue

            match                 = matches[0]
            target_motif_residues = []

            for i in range(match.start, match.end):
                target_motif_residues.append(idx_auth_map_tgt[i])

            all_motifs_residues.extend(target_motif_residues)

        predicted_chains.append(
            PredictedResiduesPolymer(
                polymer_class=nomenclature_class,
                source=PredictionSource(
                    source_seq=seq_src,
                    auth_asym_id=nbr_polymer["auth_asym_id"],
                    source_bound_residues=[
                        ResidueSummary(
                            auth_seq_id         = residue[0],
                            resname             = residue[1],
                            parent_auth_asym_id = nbr_polymer["auth_asym_id"],
                            label_seq_id        = None,
                            full_id             = None,
                        )
                        for residue in nbr_polymer["bound_residues"]
                    ],
                ),
                target=PredictionTarget(
                    target_seq            = seq_tgt,
                    auth_asym_id          = target_polymer.auth_asym_id,
                    target_bound_residues = [
                        ResidueSummary(
                            auth_seq_id         = residue.get_id()[1],
                            resname             = ResidueSummary.three_letter_code_to_one(residue.resname),
                            label_seq_id        = None,
                            parent_auth_asym_id = target_polymer.auth_asym_id,
                            full_id             = None,
                        )
                        for residue in all_motifs_residues
                    ],
                ),
            )
        )

    # ! at this point we have collected all the source polymers and corresponding target polymers.
    chains: list[BindingSiteChain] = []

    for c in predicted_chains:
        poly = RibosomeOps(target_rcsb_id).get_poly_by_auth_asym_id(
            c.target.auth_asym_id
        )
        if poly == None:
            raise ValueError("Polymer not found in target structure")
        chains.append(
            BindingSiteChain(
                **poly.model_dump(),
                bound_residues=c.target.target_bound_residues,
            )
        )

    return LigandTransposition(
        constituent_chains=predicted_chains,  # this is the info about each individual pair of polymers manipulations
        source=source_rcsb_id,
        target=target_rcsb_id,
        purported_binding_site=BindingSite(  # this is the result, the datastructure that gets sent the fronted
            chains=chains,
            ligand=binding_site.ligand,
            radius=binding_site.radius,
            source=binding_site.source,
        ),
    )


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
