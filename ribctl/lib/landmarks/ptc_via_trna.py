from pprint import pprint
from typing import List, NewType, Tuple, TypeVar
import typing
from ribctl.etl.assets_global import GlobalView
from ribctl.etl.assets_structure import StructureAssetPaths
from ribctl.ribosome_ops import RibosomeOps
import pickle
from Bio.PDB.Residue import Residue
import copy
import typing
from Bio.PDB.Residue import Residue
from Bio.PDB.Chain import Chain
import numpy as np
from ribctl.ribosome_ops import RibosomeOps
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection
from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import SequenceMappingContainer
from ribctl.lib.libtax import Taxid
from ribctl.lib.schema.types_binding_site import ResidueSummary
from scipy.spatial.distance import pdist, squareform

REFERENCE_MITO_STRUCTURE_TRNA_RRNA = ("7A5F", "24", "A3")
REFERENCE_ARCHAEA_STRUCTURE_TRNA_RRNA = ("8HKY", "APTN", "A23S")
REFERENCE_BACTERIA_STRUCTURE_TRNA_RRNA = ("8UD8", "1x", "1A")
REFERENCE_EUKARYA_STRUCTURE_TRNA_RRNA = ("8CCS", "Bb", "AA")


def find_closest_pair(points: np.ndarray):
    points = np.asarray(points)
    if len(points) < 2:
        raise ValueError("Array must contain at least 2 points")

    distances = pdist(points)
    distance_matrix = squareform(distances)
    i, j = np.triu_indices(len(points), k=1)
    min_idx = np.argmin(distances)
    point1_idx = i[min_idx]
    point2_idx = j[min_idx]

    closest_point1 = points[point1_idx]
    closest_point2 = points[point2_idx]
    min_distance = distances[min_idx]

    return closest_point1, closest_point2, min_distance


def PTC_reference_residues(
    ribosome_type: typing.Literal["euk", "bact", "arch", "mito"]
) -> tuple[list[Residue], Chain, tuple[str, str, str]]:
    match ribosome_type:
        case "euk":
            ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = (
                REFERENCE_EUKARYA_STRUCTURE_TRNA_RRNA
            )
        case "bact":
            ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = (
                REFERENCE_BACTERIA_STRUCTURE_TRNA_RRNA
            )
        case "arch":
            ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = (
                REFERENCE_ARCHAEA_STRUCTURE_TRNA_RRNA
            )
        case "mito":
            ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = (
                REFERENCE_MITO_STRUCTURE_TRNA_RRNA
            )
        case _:
            raise ValueError("Invalid ribosome type")

    print(
        "\t Seeking the LSU rRNA[{}] residues in the  vicinity of tRNA[{}] chain's C-terminus in [{}]".format(
            ref_rrna_aaid, ref_trna_aaid, ref_rcsb_id
        )
    )

    mmcif_struct = RibosomeOps(ref_rcsb_id).assets.biopython_structure()[0]

    def trna_cterm_pos() -> np.ndarray:
        trnaChain: Chain = mmcif_struct[ref_trna_aaid]
        c_terminus: Residue = list(
            filter(
                lambda x: ResidueSummary.filter_noncanonical(x.resname), [*trnaChain]
            )
        )[-1]
        return c_terminus.center_of_mass()

    rrrna = mmcif_struct[ref_rrna_aaid]
    atoms = Selection.unfold_entities(rrrna, "A")
    ns = NeighborSearch(atoms)
    nearby_residues = ns.search(trna_cterm_pos(), 10, "R")

    return (
        list(
            filter(
                lambda x: ResidueSummary.filter_noncanonical(x.resname), nearby_residues
            )
        ),
        rrrna,
        (ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid),
    )


def pickle_ref_ptc_data(ref_data: dict, output_file: str):
    try:
        with open(output_file, "wb") as f:
            pickle.dump(ref_data, f, protocol=pickle.HIGHEST_PROTOCOL)
            print("Saved {}".format(output_file))
        return True

    except Exception as e:
        print(f"Error pickling residues: {str(e)}")
        return False


def unpickle_residue_array(input_file: str) -> dict | None:
    try:
        with open(input_file, "rb") as f:
            data_dict = pickle.load(f)
        return data_dict
    except Exception as e:
        print(f"Error unpickling residues: {str(e)}")
        return None


def produce_ptc_references():
    for ribosome_type in ["mito", "euk", "arch", "bact"]:
        residues, chain, meta = PTC_reference_residues(ribosome_type)
        ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid = meta
        _ = {
            "nearest_residues": residues,
            "chain": chain,
            "ref_rcsb_id": ref_rcsb_id,
            "ref_trna_aaid": ref_trna_aaid,
            "ref_rrna_aaid": ref_rrna_aaid,
        }
        outpath = GlobalView.ptc_references(ribosome_type)
        pickle_ref_ptc_data(_, outpath)


def get_ptc_reference(ribosome_type: typing.Literal["mito", "euk", "arch", "bact"]):
    cached_name = GlobalView.ptc_references(ribosome_type)
    return unpickle_residue_array(cached_name)


def PTC_location(target_rcsb_id: str) -> Tuple[np.ndarray, list[Residue]]:
    """
    # Get PTC in @target_rcsb_id by way of mapping a reference PTC in a given mitochondrial structure
    #"""
    RO = RibosomeOps(target_rcsb_id)
    tax_id = RO.taxid
    match Taxid.superkingdom(tax_id):
        case "archaea":
            ribosome_type = "arch"
        case "bacteria":
            ribosome_type = "bact"
        case "eukaryota":
            ribosome_type = "euk"
        case _:
            raise ValueError("Invalid taxid")

    if RO.profile.mitochondrial:
        ribosome_type = "mito"

    data_dict = get_ptc_reference(ribosome_type)
    if data_dict is None:
        raise IndexError("Reference file doesn't exist. It should")

    ref_residues: list[Residue] = data_dict["nearest_residues"]
    ref_chain: Chain = data_dict["chain"]
    # data_dict['ref_rcsb_id'     ]
    # data_dict['ref_trna_aaid'   ]
    # data_dict['ref_rrna_aaid'   ]

    mmcif_struct_tgt = RO.assets.biopython_structure()[0]
    LSU_RNA_tgt_aaid = RO.get_LSU_rRNA().auth_asym_id
    LSU_RNA_tgt: Chain = mmcif_struct_tgt[LSU_RNA_tgt_aaid]
    print("got target", [LSU_RNA_tgt.child_dict])
    print("got target", LSU_RNA_tgt)

    _, _, motifs = map_motifs(
        SequenceMappingContainer(ref_chain),
        SequenceMappingContainer(LSU_RNA_tgt),
        [ResidueSummary.from_biopython_residue(r) for r in ref_residues],
        "-",
        True,
    )

    # RO     = RibosomeOps(target_rcsb_id)
    # mmcif_struct_tgt = RO.biopython_structure()[0]
    # LSU_RNA_tgt_aaid = RO.get_LSU_rRNA().auth_asym_id
    # LSU_RNA_tgt:Chain      = mmcif_struct_tgt["72"]
    # nums  = [2602,2451 ]
    # residues =LSU_RNA_tgt.child_list
    # x :Residue= residues[1]
    # motifs = list(filter( lambda x: int(x.full_id[3][1]) in nums, LSU_RNA_tgt.child_list ))
    # pprint(motifs)

    (p1, p2, dist) = find_closest_pair([r.center_of_mass() for r in motifs])
    return (p1 + p2) / 2, motifs
