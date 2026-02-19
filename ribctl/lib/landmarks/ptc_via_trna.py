# ribctl/lib/landmarks/ptc_via_trna.py

from __future__ import annotations

import typing
import pickle

import numpy as np
from Bio.PDB import Selection
from Bio.PDB.Chain import Chain
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from scipy.spatial.distance import pdist, squareform

from ribctl.global_ops import GlobalOps
from ribctl.lib.exceptions import SkipAsset
from ribctl.lib.libbsite import map_motifs
from ribctl.lib.libseq import SequenceMappingContainer
from ribctl.lib.libtax import Taxid
from ribctl.lib.schema.types_binding_site import ResidueSummary
from ribctl.lib.schema.types_ribosome import PTCInfo
from ribctl.ribosome_ops import RibosomeOps


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
    ribosome_type: typing.Literal["euk", "bact", "arch", "mito"],
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
        canon = list(
            filter(lambda x: ResidueSummary.is_canonical(x.resname), [*trnaChain])
        )
        if not canon:
            raise SkipAsset(
                f"Reference tRNA chain {ref_trna_aaid} in {ref_rcsb_id} contains no canonical residues"
            )
        c_terminus: Residue = canon[-1]
        return c_terminus.center_of_mass()

    rrrna = mmcif_struct[ref_rrna_aaid]
    atoms = Selection.unfold_entities(rrrna, "A")
    ns = NeighborSearch(atoms)
    nearby_residues = ns.search(trna_cterm_pos(), 10, "R")

    filtered = list(
        filter(lambda x: ResidueSummary.is_canonical(x.resname), nearby_residues)
    )
    if not filtered:
        raise SkipAsset(
            f"No canonical nearby rRNA residues found in reference {ref_rcsb_id} ({ribosome_type})"
        )

    return (filtered, rrrna, (ref_rcsb_id, ref_trna_aaid, ref_rrna_aaid))


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
        outpath = GlobalOps.ptc_references(ribosome_type)
        pickle_ref_ptc_data(_, outpath)


def get_ptc_reference(ribosome_type: typing.Literal["mito", "euk", "arch", "bact"]):
    cached_name = GlobalOps.ptc_references(ribosome_type)
    return unpickle_residue_array(cached_name)


def _infer_ribosome_type(
    ro: RibosomeOps,
) -> typing.Literal["mito", "euk", "arch", "bact"]:
    """
    Determine which reference bucket to use.
    """
    if ro.profile.mitochondrial:
        return "mito"

    tax_id = ro.taxid
    match Taxid.superkingdom(tax_id):
        case "archaea":
            return "arch"
        case "bacteria":
            return "bact"
        case "eukaryota":
            return "euk"
        case _:
            raise SkipAsset(f"Unsupported/unknown superkingdom for taxid={tax_id}")


def _has_lsu(profile) -> bool:
    """
    Prefer the explicit subunit_presence annotation when available.
    """
    try:
        if profile.subunit_presence:
            return "lsu" in profile.subunit_presence
    except Exception:
        pass
    # If absent/unreliable, don't hard-fail here; we'll attempt LSU rRNA lookup next.
    return True


def _try_get_lsu_rrna_poly(ro: RibosomeOps):
    """
    Try to get LSU rRNA across all assembly_ids present in profile.rnas.
    This avoids hard-coding assembly=0 assumptions.
    """
    profile = ro.profile
    assembly_ids = set()
    for rna in getattr(profile, "rnas", None) or []:
        try:
            assembly_ids.add(int(rna.assembly_id))
        except Exception:
            pass
    # Always try 0 first as a common case
    ordered = [0] + sorted(a for a in assembly_ids if a != 0)

    last_err = None
    for aid in ordered:
        try:
            return ro.get_LSU_rRNA(assembly=aid)
        except Exception as e:
            last_err = e
            continue

    raise SkipAsset(
        f"No LSU rRNA found in structure (tried assemblies {ordered}): {last_err}"
    )


def PTC_location(target_rcsb_id: str) -> PTCInfo:
    """
    Get PTC in @target_rcsb_id by mapping reference PTC-adjacent residues
    from a reference rRNA chain onto the target LSU rRNA chain.
    """
    RO = RibosomeOps(target_rcsb_id)
    profile = RO.profile

    # Fast skip on SSU-only if annotation exists
    if not _has_lsu(profile):
        raise SkipAsset("No LSU annotated in subunit_presence; cannot compute PTC")

    ribosome_type = _infer_ribosome_type(RO)

    data_dict = get_ptc_reference(ribosome_type)
    if data_dict is None:
        raise RuntimeError(
            f"PTC reference cache missing/unreadable for ribosome_type={ribosome_type}. "
            f"Expected file at {GlobalOps.ptc_references(ribosome_type)}"
        )

    ref_residues: list[Residue] = data_dict["nearest_residues"]
    ref_chain: Chain = data_dict["chain"]

    mmcif_struct_tgt = RO.assets.biopython_structure()[0]

    # Locate LSU rRNA polymer and chain robustly across assemblies
    lsu_poly = _try_get_lsu_rrna_poly(RO)
    LSU_RNA_tgt_aaid = lsu_poly.auth_asym_id

    try:
        LSU_RNA_tgt: Chain = mmcif_struct_tgt[LSU_RNA_tgt_aaid]
    except KeyError:
        # Some files place chains in different models; scan models
        found = None
        for model in RO.assets.biopython_structure():
            try:
                if LSU_RNA_tgt_aaid in model.child_dict:
                    found = model[LSU_RNA_tgt_aaid]
                    break
            except Exception:
                continue
        if found is None:
            raise SkipAsset(
                f"LSU rRNA chain {LSU_RNA_tgt_aaid} not found in mmCIF models"
            )
        LSU_RNA_tgt = found

    # Map reference residues to target residues via sequence mapping
    _, _, motifs = map_motifs(
        SequenceMappingContainer(ref_chain),
        SequenceMappingContainer(LSU_RNA_tgt),
        [ResidueSummary.from_biopython_residue(r) for r in ref_residues],
        "-",
        False,
    )

    if motifs is None or len(motifs) == 0:
        raise SkipAsset(
            f"PTC motif mapping yielded 0 residues for {target_rcsb_id} (type={ribosome_type})"
        )

    # Compute center of mapped residues
    coords = []
    for r in motifs:
        try:
            coords.append(r.center_of_mass())
        except Exception:
            continue

    if len(coords) == 0:
        raise SkipAsset(
            f"Mapped motifs present but no valid COM coordinates for {target_rcsb_id}"
        )

    if len(coords) == 1:
        center = coords[0]
    else:
        (p1, p2, dist) = find_closest_pair(coords)
        center = (p1 + p2) / 2

    return PTCInfo(
        location=center.tolist(),
        residues=list(map(ResidueSummary.from_biopython_residue, motifs)),
    )
