import json
import os
import operator
import argparse
import itertools
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from ribctl.lib.types.types_ribosome import Polymer, RibosomeStructure
from ribctl.lib.types.ligands.types_binding_site import AMINO_ACIDS, NUCLEOTIDES, BindingSite, BindingSiteChain, ResidueSummary
from ribctl.lib import RIBETL_DATA
from ribctl.lib import utils


def get_polymer_residues(auth_asym_id: str, struct: Structure) -> list[Residue]:
    c: Chain = struct[0][auth_asym_id]
    return [*c.get_residues()]

def get_polymer_nbrs(
      residues      : list[Residue],
      struct        : Structure,
    ) -> BindingSite: 

    """KDTree search the neighbors of a given list of residues(which constitue a ligand)
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid         = struct.get_id().upper()
    ns            = NeighborSearch(list(struct.get_atoms()))
    nbr_residues  = []
    parent_strand = residues[0].get_parent().id if len(residues) > 0 else ...

    for poly_res in residues:
        for atom in poly_res.child_list:
            nbr_residues.extend(ns.search(atom.get_coord(), 10, level='R'))

    nbr_residues = list(set([* map(ResidueSummary.from_biopython_residue, nbr_residues)]))
    nbr_residues = list(filter(lambda res: res.parent_auth_asym_id != parent_strand and res.get_resname() in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _: _.parent_auth_asym_id, nbr_residues)))

    

    profile       = RibosomeStructure.from_json_profile(utils.open_structure(pdbid, 'json'))
    poly_entities:list[Polymer] = [*map(lambda _: Polymer(**_.dict()), itertools.chain(profile.proteins, profile.rnas if profile.rnas else []))]

    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity.auth_asym_id:

                    nbr_dict[c] = BindingSiteChain(
                        **poly_entity.dict(),
                        residues     = sorted([residue for residue in nbr_residues if residue.parent_auth_asym_id == c], key=operator.attrgetter('seqid')))

    return BindingSite(nbr_chains=nbr_dict)

def get_ligand_nbrs(
      ligand_residues: list[Residue],
      struct         : Structure,
    ) -> BindingSite : 
    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid = struct.get_id().upper()

    profile = RibosomeStructure.from_json_profile(utils.open_structure(pdbid, 'json'))
    poly_entities = itertools.chain(profile.proteins, profile.rnas if profile.rnas else [])

    pdbid        = struct.get_id()
    ns           = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
            nbr_residues.extend(ns.search(atom.get_coord(), 10, level='R'))

    nbr_residues = list(set([* map(ResidueSummary.from_biopython_residue, nbr_residues)]))

    # Filter the ligand itself, water and other special residues
    nbr_residues = list(filter(lambda resl: resl.resname in [*AMINO_ACIDS.keys(),  *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.get_parent_auth_asym_id(), nbr_residues)))
    

    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity.auth_asym_id:
                    nbr_dict[c] = BindingSiteChain(
                        **poly_entity.dict(),
                        residues=sorted(
                            [residue for residue in nbr_residues if residue.get_parent_auth_asym_id() == c],
                            key=operator.attrgetter('seqid')
                        ))

    return BindingSite(nbr_chains=nbr_dict)

def get_ligand_residue_ids(ligchemid: str, struct: Structure) -> list[Residue]:
    ligandResidues: list[Residue] = list(filter(lambda x: x.get_resname() == ligchemid, list(struct.get_residues())))
    return ligandResidues

def struct_liglike_ids(struct_profile:RibosomeStructure) -> list[Polymer]:
    """Given an rcsb id, open the profile for the corresponding structure
    and return references to all polymers marked ligand-like"""
    polymers =  itertools.chain(struct_profile.rnas if struct_profile.rnas != None else [] , struct_profile.proteins)
    return list(filter(lambda poly : poly.ligand_like == True, polymers))

def struct_ligand_ids(pdbid: str, profile:RibosomeStructure) -> list[tuple]:
    pdbid = pdbid.upper()
    if not profile.ligands: return []
    _ = [* map(lambda x: (x.chemicalId, x.chemicalName),profile.ligands)] 
    return [ ] if len(_) < 1 else [* filter(lambda k: "ion" not in k[1].lower(), _)] 

def render_liglike_polymer(rcsb_id:str, auth_asym_id:str, structure:Structure, WRITE:bool=False)->BindingSite:
    residues: list[Residue] = get_polymer_residues(auth_asym_id, structure)
    binding_site_polymer: BindingSite = get_polymer_nbrs(residues, structure )

    outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'POLYMER_{auth_asym_id}.json')
    if WRITE:
        with open(outfile_json, 'w') as outfile:
            json.dump(binding_site_polymer.dict(), outfile, indent=4)
            print("Wrote: ", outfile_json)
    else:
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)

    return binding_site_polymer

def render_ligand(rcsb_id:str,chemicalId:str, structure:Structure, WRITE:bool=False)->BindingSite:
    chemicalId = chemicalId.upper()
    residues: list[Residue] = get_ligand_residue_ids(chemicalId, structure)
    binding_site_ligand: BindingSite   = get_ligand_nbrs(residues, structure)

    outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'LIGAND_{chemicalId}.json')
    if WRITE:
        with open(outfile_json, 'w') as outfile:
            json.dump(binding_site_ligand.dict(), outfile, indent=4)
            print("Wrote: ", outfile_json)

    elif (os.path.isfile(outfile_json)):
        print("Exists already: ", outfile_json)

    return binding_site_ligand

if __name__ == "__main__":

    parser = argparse. ArgumentParser(description='Split structure into constituent polymers and inject new nomencalture into the .cif file')
    parser.add_argument ('-s'     , '--structure', type=   str   , required=True                                                            )
    parser.add_argument ('--save' ,action          ='store_true'                                                    )
    
    args  = parser.parse_args()
    PDBID = args.structure.upper()

    _structure_cif_handle :Structure = utils.open_structure(PDBID,'cif')
    struct_profile_handle       = RibosomeStructure.parse_obj(utils.open_structure(PDBID,'json'))

    liglike_polys = struct_liglike_ids(struct_profile_handle)
    ligands       = struct_ligand_ids(PDBID, struct_profile_handle)


    for polyref in liglike_polys:
        render_liglike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, args.save)

    for l in ligands:
        render_ligand(PDBID, l[0], _structure_cif_handle, args.save)
