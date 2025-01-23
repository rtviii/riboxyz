from pprint import pprint
import sys
from ribctl.lib.libseq import SequenceMappingContainer, SeqPairwise, map_motifs
sys.path.append("/home/rtviii/dev/riboxyz")
# from neo4j_ribosome.db_lib_reader import Neo4jReader
import operator
import json
import warnings
from Bio import ( BiopythonDeprecationWarning, )
import pandas as pd
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonDeprecationWarning)
from ribctl.lib.schema.types_binding_site import (
    BindingSiteChain,
    LigandTransposition,
    ResiduesMapping,
    PredictionSource,
    PredictionTarget,
)

from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.schema.types_binding_site import (
    BindingSite,
    BindingSiteChain,
    ResidueSummary,
)
import os
from collections import defaultdict
from ribctl import ASSETS_PATH

def lig_get_chemical_categories():
    ligands_classification_path = os.path.join(
        ASSETS_PATH, "ligands", "ligand_chemical_categories.csv"
    )
    df = pd.read_csv(ligands_classification_path)
    reverse_dict = defaultdict(list)
    for _, row in df.iterrows():
        ligand = row["Ligand"]
        category = row["Category"]
        reverse_dict[category].append(ligand)

    return dict(reverse_dict)

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
    assemblies = ( RibosomeOps(struct.get_id().upper()).profile.get_polymers_by_assembly() )
    # If there are two or more assemblies, delete chains belonging to all but the first one.
    if len(assemblies.items()) > 1:
        for i in range(len(assemblies.items()) - 1):
            for chain_aaid in [*assemblies.items()][i + 1][1]:
                md.detach_child(chain_aaid)

    ns              = NeighborSearch(list(struct.get_atoms()))

    nbr_residues    = []
    ligand_residues = list(filter(lambda x: x.get_resname() == lig_chemid, list(struct.get_residues())) )

    #** So, for the case where there are multiple ligands, we only take the first one.
    #** Otherwise there is possibility of "disjoint" pockets being detected which is not convenient.
    #** The identification system could be better at the structural level on my or mmcif side.
    for lig_res in ligand_residues[:1]:
        for atom in lig_res.child_list:
            found_nbrs = ns.search(atom.get_coord(), radius, level="R")
            found_nbrs = [r for r in found_nbrs if r.get_resname() != lig_chemid]
            nbr_residues.extend(found_nbrs)

    nbr_residues: list[Residue] = list(set(nbr_residues))
    nbr_chains = []

    nbr_residues_by_chain_aaid = {}

    for residue in nbr_residues:

        parent_chain = residue.get_parent()
        auth_asym_id = parent_chain.get_id()

        if auth_asym_id not in nbr_residues_by_chain_aaid:
            nbr_residues_by_chain_aaid[auth_asym_id] = [residue]
        else:
            nbr_residues_by_chain_aaid[auth_asym_id].append(residue)


    RO = RibosomeOps(struct.get_id().upper())

    pprint(list(nbr_residues_by_chain_aaid.keys()))
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
    _structure_cif_handle = RibosomeOps(rcsb_id).assets.biopython_structure()
    binding_site_ligand = get_lig_bsite(chemicalId, _structure_cif_handle, radius)

    if save:
        with open(RibosomeOps(rcsb_id).assets.paths.binding_site(chemicalId), "w") as f:
            json.dump(binding_site_ligand.model_dump(), f)

    return binding_site_ligand

def bsite_transpose(
    source_rcsb_id: str,
    target_rcsb_id: str,
    binding_site  : BindingSite,
    save          : bool = False,
    verbose       : bool=False

) -> LigandTransposition:

    source_rcsb_id, target_rcsb_id = source_rcsb_id.upper(), target_rcsb_id.upper()
    target_struct = RibosomeOps(target_rcsb_id).assets.biopython_structure()
    source_struct = RibosomeOps(source_rcsb_id).assets.biopython_structure()

    source_ops = RibosomeOps(source_rcsb_id)
    # source_profile = source_ops.profile()

    target_ops = RibosomeOps(target_rcsb_id)
    # target_profile = target_ops.profile()

    chain_mappings      = []
    bindign_site_chains = []

    for source_polymer in binding_site.chains:

        source_polymer = BindingSiteChain.model_validate(source_polymer)

        #! Skip if no nomenclature present ( can't do anything with it )
        if len(source_polymer.nomenclature) < 1:
            continue

        target_polymer = target_ops.get_poly_by_polyclass( source_polymer.nomenclature[0], 0 )
        if target_polymer == None:
            continue

        bpchain_source = SequenceMappingContainer(source_struct[0][source_polymer.auth_asym_id])
        bpchain_target = SequenceMappingContainer(target_struct[0][target_polymer.auth_asym_id])

        primary_seq_source, primary_seq_target, tgt_bound_residues =  map_motifs(bpchain_source, bpchain_target, source_polymer.bound_residues, source_polymer.nomenclature[0], verbose)

        polymer_pair = ResiduesMapping(
            polymer_class = source_polymer.nomenclature[0],
            source        = PredictionSource(
                source_seq            = primary_seq_source,
                auth_asym_id          = source_polymer.auth_asym_id,
                source_bound_residues = [
                    ResidueSummary(
                        auth_seq_id   = residue.auth_seq_id,
                        label_comp_id = residue.label_comp_id,
                        auth_asym_id  = source_polymer.auth_asym_id,
                        label_seq_id  = None,
                        full_id       = None,
                        rcsb_id       = source_rcsb_id,
                    )
                    for residue in source_polymer.bound_residues
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
                    ) for residue in tgt_bound_residues
                ],
            ),
        )
        chain_mappings.append(polymer_pair)
        bindign_site_chains.append(
            BindingSiteChain(
                **target_polymer.model_dump(),
                bound_residues=[
                    ResidueSummary(
                        auth_seq_id=residue.get_id()[1],
                        label_comp_id=residue.resname,
                        label_seq_id=None,
                        auth_asym_id=target_polymer.auth_asym_id,
                        full_id=None,
                        rcsb_id=target_rcsb_id,
                    )
                    for residue in tgt_bound_residues
                ],
            )
        )

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


