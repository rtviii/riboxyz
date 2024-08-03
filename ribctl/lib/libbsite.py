import json
import operator
import argparse
import itertools
from pprint import pprint
from typing import Optional
import re
from typing import List
import warnings
from Bio import (
    BiopythonExperimentalWarning,
    BiopythonWarning,
    BiopythonDeprecationWarning,
)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import pairwise2
from ribctl.lib.schema.types_binding_site import (
    BindingSiteChain,
    LigandTransposition,
    PredictedResiduesPolymer,
)
from ribctl.lib.schema.types_ribosome import (
    Polymer,
    PolymerClass,
    PolynucleotideClass,
    RibosomeStructure,
)
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from ribctl.etl.etl_assets_ops import Assets, RibosomeOps
from ribctl.lib.schema.types_ribosome import (
    NonpolymericLigand,
    Polymer,
    RibosomeStructure,
)
from ribctl.lib.schema.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
    BindingSite,
    BindingSiteChain,
    ResidueSummary,
)
from ribctl.lib import utils


def get_polymer_residues(auth_asym_id: str, struct: Structure) -> list[Residue]:
    c: Chain = struct[0][auth_asym_id]
    return [*c.get_residues()]


def get_lig_bsite(
    ligand_residues: list[Residue],
    struct: Structure,
    radius: Optional[float] = 5,

) -> BindingSite:
    """KDTree search the neighbors of a given list of residues(which constitue a ligand)
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues.
    """

    RO =RibosomeOps(struct.get_id().upper())
    profile       = RO.profile()
    ns            = NeighborSearch(list(struct.get_atoms()))
    nbr_residues  = []

    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
            nbr_residues.extend(ns.search(atom.get_coord(), radius, level="R"))

    nbr_residues = list( set([*map(ResidueSummary.from_biopython_residue, nbr_residues)]) )

    # Filter the ligand itself, water and other special residues
    nbr_residues = list(
        filter(
            lambda resl: resl.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES],
            nbr_residues,
        )
    )

    nbr_chains = []
    chain_auth_asym_ids = list(set(map(lambda _: _.get_parent_auth_asym_id(), nbr_residues)))

    for c in chain_auth_asym_ids:

            polymer= RO.get_poly_by_auth_asym_id(c)
            if polymer == None:
                raise ValueError(f"Polymer with auth_asym_id {c} not found in structure")
            nbr_chains.append(BindingSiteChain(
                **polymer.model_dump(),
                residues=sorted(
                    [
                        residue
                        for residue in nbr_residues
                        if residue.get_parent_auth_asym_id() == c
                    ],
                    key=operator.attrgetter("seqid"),
                ),
            ))
    return BindingSite.model_validate({"chains":nbr_chains })


def get_ligand_residue_ids(ligchemid: str, struct: Structure) -> list[Residue]:
    ligandResidues: list[Residue] = list(
        filter(lambda x: x.get_resname() == ligchemid, list(struct.get_residues()))
    )
    return ligandResidues


def struct_ligand_ids(
    pdbid: str, profile: RibosomeStructure
) -> list[NonpolymericLigand]:
    """
    we identify ligands worth interest by their having drugbank annotations (ions and water are excluded)
    """

    pdbid = pdbid.upper()
    ligs = []

    if not profile.nonpolymeric_ligands:
        return []
    for lig in profile.nonpolymeric_ligands:
        if lig.nonpolymer_comp.drugbank == None:
            continue
        else:
            ligs.append(lig)
    return ligs


def bsite_extrarbx_polymer(auth_asym_id: str, structure: Structure) -> BindingSite:

    residues: list[Residue] = get_polymer_residues(auth_asym_id, structure)
    binding_site_polymer: BindingSite = get_polymer_nbrs(residues, structure)
    return binding_site_polymer


def bsite_ligand(
    chemicalId: str, rcsb_id: str, radius: Optional[float] = 5, save: bool = False
) -> BindingSite:
    chemicalId = chemicalId.upper()
    _structure_cif_handle = RibosomeOps(rcsb_id).biopython_structure()
    residues: list[Residue] = get_ligand_residue_ids(chemicalId, _structure_cif_handle)
    binding_site_ligand: BindingSite = get_lig_bsite(
        residues, _structure_cif_handle, radius
    )

    if save:
        with open(RibosomeOps(rcsb_id).paths.binding_site(chemicalId), "w") as f:
            json.dump(binding_site_ligand.model_dump(), f)

    return binding_site_ligand


#! Transposition methods
class SeqMatch:
    def __init__(self, sourceseq: str, targetseq: str, source_residues: list[int]):
        """A container for origin and target sequences when matching the resiudes of a ligand binding site
        to another protein's sequence through BioSeq's Align
        """

        # * Computed indices of the ligand-facing in the source sequence.
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


# TODO: This signature has to contain info for both soruce and target rcsb_id. This is retarded as it is rn.
def init_transpose_ligand(
    source_struct: str,
    target_struct: str,
    binding_site: BindingSite,
    save: bool = False,
) -> LigandTransposition:

    by_polymer_class_source_polymers = {}
    for nbr_polymer in binding_site:
        nbr_polymer = BindingSiteChain.model_validate(nbr_polymer)
        if len(nbr_polymer.nomenclature) < 1:
            continue
        else:
            print(
                "auth asym id "
                + nbr_polymer.auth_asym_id
                + " has  nomenclature: "
                + nbr_polymer.nomenclature[0]
            )
            by_polymer_class_source_polymers[nbr_polymer.nomenclature[0].value] = {
                "seq": nbr_polymer.entity_poly_seq_one_letter_code_can,
                "auth_asym_id": nbr_polymer.auth_asym_id,
                "ids": [
                    resid for resid in [*map(lambda x: x.seqid, nbr_polymer.residues)]
                ],
            }

    # ! at this point we have collected all the source polymers, their sequences and residue ids participating in the binding site
    by_class_target_polymers = {}
    for nomenclature_class, nbr_polymer in by_polymer_class_source_polymers.items():
        target_polymer = RibosomeOps(target_struct).get_poly_by_polyclass(
            nomenclature_class
        )
        if target_polymer == None:
            continue
        tgt_poly_seq = target_polymer.entity_poly_seq_one_letter_code_can
        tgt_poly_auth_asym_id = target_polymer.auth_asym_id
        by_class_target_polymers[nomenclature_class] = {
            "seq": tgt_poly_seq,
            "auth_asym_id": tgt_poly_auth_asym_id,
        }

    predicted_chains: list[PredictedResiduesPolymer] = []

    for nomenclature_class, seqstats in by_polymer_class_source_polymers.items():
        if nomenclature_class not in by_class_target_polymers:
            continue

        src_ids = by_polymer_class_source_polymers[nomenclature_class]["ids"]
        src = by_polymer_class_source_polymers[nomenclature_class]["seq"]
        tgt = by_class_target_polymers[nomenclature_class]["seq"]

        sq = SeqMatch(src, tgt, src_ids)

        src_aln = (
            sq.src_aln
        )  # <--- aligned source      sequence (with                        gaps)
        tgt_aln = (
            sq.tgt_aln
        )  # <--- aligned tgt         sequence (with                        gaps)
        aln_ids = (
            sq.aligned_ids
        )  # <--- ids     corrected   for                                   gaps
        tgt_ids = (
            sq.tgt_ids
        )  # <--- ids     backtracted to the target polymer (accounting for gaps)
        
        predicted_chains.append(PredictedResiduesPolymer.model_validate(
            {
                "polymer_class": nomenclature_class,
                "source": {
                    "source_seq": src,
                    "source_seq_ids": src_ids,
                    "auth_asym_id": by_polymer_class_source_polymers[
                        nomenclature_class
                    ].auth_asym_id,
                },
                "target": {
                    "target_seq": tgt,
                    "target_seq_ids": tgt_ids,
                    "auth_asym_id": by_class_target_polymers[nomenclature_class][
                        "auth_asym_id"
                    ],
                },
                "alignment": {
                    "aligned_ids": aln_ids,
                    "source_seq_aligned": src_aln,
                    "target_seq_aligned": tgt_aln,
                },
            }
        ))

    purported_binding_site: BindingSite = []

    for chain in predicted_chains:
        poly = RibosomeOps(target_struct).get_poly_by_auth_asym_id(chain.target.auth_asym_id) 
        if poly == None: raise ValueError("Polymer not found in target structure")

        purported_binding_site.append(BindingSiteChain(
            **poly.model_dump(),
            residues=[
                ResidueSummary(
					seqid               = seqid,
					resname             = chain.target.target_seq[seqid],
					parent_auth_asym_id = chain.target.auth_asym_id,
					full_id             = None
                ) for  seqid in chain.target.target_seq_ids
            ],
        ))

    return LigandTransposition(
        constituent_chains     = predicted_chains,
        source                 = source_struct,
        target                 = target_struct,
        purported_binding_site = purported_binding_site,
    )
