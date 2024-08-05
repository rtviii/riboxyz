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
    BiopythonDeprecationWarning,
)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
    from Bio import pairwise2
from ribctl.lib.schema.types_binding_site import (
    BindingSiteChain,
    LigandTransposition,
    PredictedResiduesPolymer,
    PredictionAlignments,
    PredictionSource,
    PredictionTarget,
)
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from ribctl.etl.etl_assets_ops import Assets, RibosomeOps
from ribctl.lib.schema.types_binding_site import (
    AMINO_ACIDS,
    NUCLEOTIDES,
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


    c:Chain = struct[0]['DA']
    print([*c.get_residues()])
    print(len([*c.get_residues()]))

    exit()
    ns = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    ligand_residues: list[Residue] = list( filter(lambda x: x.get_resname() == lig_chemid, list(struct.get_residues())) )


    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
            found_nbrs = ns.search(atom.get_coord(), radius, level="R")
            nbr_residues.extend(found_nbrs)

    nbr_residues:list[Residue] = list(set(nbr_residues))


    for residue in nbr_residues:
        print(residue.get_id(), "\t", residue.get_segid(), " ]", residue.get_full_id(), residue.resseq)

    pprint(dir( nbr_residues[0] ))
    exit()
    nbr_residues = list( set([*map(ResidueSummary.from_biopython_residue, nbr_residues)]) )


    pprint("Searchin neighbors withing radius of {}".format(radius))
    pprint("??????????????????????????????????????")
    exit()

    # Filter the ligand itself, water and other special residues
    nbr_residues = list( filter( lambda resl: resl.resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues, ) )





    pprint(nbr_residues)
    pprint("??????????????????????????????????????")
    exit()

    nbr_chains = []
    chain_auth_asym_ids = list( set(map(lambda _: _.get_parent_auth_asym_id(), nbr_residues)) )

    RO = RibosomeOps(struct.get_id().upper())
    for c in chain_auth_asym_ids:
        polymer = RO.get_poly_by_auth_asym_id(c)

        if polymer == None:
            raise ValueError(
                f"Polymer with auth_asym_id {c} not found in structure. Logic error."
            )
        bound_residues = sorted(
            [
                residue
                for residue in nbr_residues
                if residue.get_parent_auth_asym_id() == c
            ],
            key=operator.attrgetter("seqid"),
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
    source_struct: str,
    target_struct: str,
    binding_site: BindingSite,
    save: bool = False,
) -> LigandTransposition:

    source_struct, target_struct = source_struct.upper(), target_struct.upper()
    source_polymers_by_poly_class = {}


    for nbr_polymer in binding_site.chains:
        nbr_polymer = BindingSiteChain.model_validate(nbr_polymer)
        if len(nbr_polymer.nomenclature) < 1:
            continue
        else:
            source_polymers_by_poly_class[nbr_polymer.nomenclature[0].value] = {
                "seq"         : nbr_polymer.entity_poly_seq_one_letter_code_can,
                "auth_asym_id": nbr_polymer.auth_asym_id,
                "ids"         : [ resid for resid in [*map(lambda x: x.seqid, nbr_polymer.bound_residues)] ],
            }

    target_polymers_by_poly_class = {}
    for nomenclature_class, nbr_polymer in source_polymers_by_poly_class.items():
        target_polymer = RibosomeOps(target_struct).get_poly_by_polyclass(nomenclature_class )
        # print("Retrieved {} from target structure".format(nomenclature_class))
        if target_polymer == None:
            continue
        tgt_poly_seq                                      = target_polymer.entity_poly_seq_one_letter_code_can
        tgt_poly_auth_asym_id                             = target_polymer.auth_asym_id
        target_polymers_by_poly_class[nomenclature_class] = {
            "seq": tgt_poly_seq,
            "auth_asym_id": tgt_poly_auth_asym_id,
        }

    # ! at this point we have collected all the source polymers and corresponding target polymers.

    predicted_chains: list[PredictedResiduesPolymer] = []

    for nomenclature_class, _ in source_polymers_by_poly_class.items():
        if nomenclature_class not in target_polymers_by_poly_class:
            continue

        # print("source_polymers_by_poly_class[nomenclature_class] yields ", source_polymers_by_poly_class[nomenclature_class])
        # pprint(src_ids)
        # exit()
        src_ids = source_polymers_by_poly_class[nomenclature_class]["ids"]
        src     = source_polymers_by_poly_class[nomenclature_class]["seq"]

        tgt = target_polymers_by_poly_class[nomenclature_class]["seq"]


        src_auth_asym_id = source_polymers_by_poly_class[nomenclature_class][ "auth_asym_id" ]
        tgt_auth_asym_id = target_polymers_by_poly_class[nomenclature_class]["auth_asym_id"]

        sq = SeqMatch(src, tgt, src_ids)

        src_aln = ( sq.src_aln )  # <--- aligned source      sequence (with                        gaps)
        tgt_aln = ( sq.tgt_aln )  # <--- aligned tgt         sequence (with                        gaps)

        aln_ids = ( sq.aligned_ids )  # <--- ids     corrected   for                                   gaps
        tgt_ids = ( sq.tgt_ids )  # <--- ids     backtracted to the target polymer (accounting for gaps)

        print("got src:", src)
        print("got src len:", len(src))
        print("got srcdis:", src_ids)
        # exit()

        predicted_chains.append(
            PredictedResiduesPolymer(
                polymer_class = nomenclature_class,
                source        = PredictionSource(
                    source_seq            = src,
                    auth_asym_id          = src_auth_asym_id,
                    source_bound_residues = [
                        ResidueSummary(
                            seqid               = src_seqid,
                            resname             = src[src_seqid],
                            parent_auth_asym_id = src_auth_asym_id,
                            full_id             = None,
                        ) for src_seqid in src_ids
                    ],
                ),
                target=PredictionTarget(
                    target_seq            = tgt,
                    auth_asym_id          = tgt_auth_asym_id,
                    target_bound_residues = [
                        ResidueSummary(
                            seqid               = tgt_seqid,
                            resname             = tgt[tgt_seqid],
                            parent_auth_asym_id = tgt_auth_asym_id,
                            full_id             = None,
                        ) for tgt_seqid in tgt_ids
                    ],
                ),
                alignment=PredictionAlignments(
                    aligned_ids        = aln_ids,
                    source_seq_aligned = src_aln,
                    target_seq_aligned = tgt_aln,
                ),
            )
            # {
            #     "source_seq"     : src,
            #     "auth_asym_id"   : by_polymer_class_source_polymers[ nomenclature_class ]["auth_asym_id"],
            #     "source_bound_residues": [
            #         ResidueSummary(
            #             seqid               = seqid,
            #             resname             = src[seqid],
            #             parent_auth_asym_id = src_auth_asym_id,
            #             full_id             = None,
            #         ) for seqid in src_ids ],
            # },
            # {
            #     "target_seq"     : tgt,
            #     "auth_asym_id"   : by_class_target_polymers[nomenclature_class][ "auth_asym_id" ],
            #     "target_bound_residues": [
            #         ResidueSummary(
            #             seqid=seqid,
            #             resname=tgt[seqid],
            #             parent_auth_asym_id=tgt_auth_asym_id,
            #             full_id=None,
            #         ) for seqid in tgt_ids ],
            # },
            # {
            #     "aligned_ids"       : aln_ids,
            #     "source_seq_aligned": src_aln,
            #     "target_seq_aligned": tgt_aln,
            # },
        )

    chains: list[BindingSiteChain] = []

    for chain in predicted_chains:
        poly = RibosomeOps(target_struct).get_poly_by_auth_asym_id(
            chain.target.auth_asym_id
        )
        if poly == None:
            raise ValueError("Polymer not found in target structure")

        chains.append(
            BindingSiteChain(
                **poly.model_dump(),
                bound_residues=chain.target.target_bound_residues,
            )
        )

    return LigandTransposition(
        constituent_chains=predicted_chains,
        source=source_struct,
        target=target_struct,
        purported_binding_site=BindingSite(
            chains=chains,
            ligand=binding_site.ligand,
            radius=binding_site.radius,
            source=binding_site.source,
        ),
    )
