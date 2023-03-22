import dataclasses
import json
import os
import operator
import argparse
import itertools
from pprint import pprint
from typing import Any
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from dataclasses import dataclass, field
from ribctl.lib.types.types_ribosome import Polymer, RibosomeStructure
from ribctl.lib import RIBETL_DATA
from ribctl.lib import utils
from ribctl.lib.types.ligands.types_binding_site import AMINO_ACIDS, NUCLEOTIDES, BindingSite, BindingSiteChain, ResidueSummary

flatten = itertools.chain.from_iterable


def __get_polymer_residues(auth_asym_id: str, struct: Structure) -> list[Residue]:
    c: Chain = struct[0][auth_asym_id]
    return [*c.get_residues()]

def __get_poly_nbrs(
      residues      : list[Residue],
      struct        : Structure,
      auth_asym_id  : str
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

    

    profile       = utils.open_structure(pdbid, 'json')
    poly_entities = [*profile['proteins'], *profile['rnas']]

    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity['auth_asym_id']:
                    nomenclature = poly_entity['nomenclature']
                    asym_id      = poly_entity['asym_ids']
                    auth_asym_id = poly_entity['auth_asym_id']
                    seq          = poly_entity['entity_poly_seq_one_letter_code']

                    nbr_dict[c] = BindingSiteChain(
                        sequence     = seq,
                        nomenclature = nomenclature,
                        asym_ids     = asym_id,
                        auth_asym_id = auth_asym_id,
                        residues     = sorted([residue for residue in nbr_residues if residue.parent_auth_asym_id == c], key=operator.attrgetter('seqid')))


    return BindingSite(__root__=nbr_dict)

def __get_ligand_nbrs(
      ligand_residues: list[Residue],
      struct         : Structure,
    ) -> BindingSite : 
    """KDTree search the neighbors of a given list of residues(which constitue a ligand) 
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid = struct.get_id().upper()

    with open(os.path.join(RIBETL_DATA, pdbid, f"{pdbid}.json"), 'rb') as strfile:
        profile = json.load(strfile)
        poly_entities = [*profile['proteins'], *profile['rnas']]

    pdbid        = struct.get_id()
    ns           = NeighborSearch(list(struct.get_atoms()))
    nbr_residues = []

    for lig_res in ligand_residues:
        for atom in lig_res.child_list:
            nbr_residues.extend(ns.search(atom.get_coord(), 10, level='R'))

    # ? Filtering phase
    # Convert residues to the dataclass (hashable), filter non-unique
    print("residues look thus:",nbr_residues)

    for n in nbr_residues:
        print(" id :", n.get_id())
        print(" id :", n.resname)
    nbr_residues = list(set([* map(ResidueSummary.from_biopython_residue, nbr_residues)]))
    pprint(nbr_residues)
    

    # Filter the ligand itself, water and other special residues
    nbr_residues = list(filter(lambda resl: resl.resname in [*AMINO_ACIDS.keys(),  *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.get_parent_auth_asym_id(), nbr_residues)))
    

    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity['auth_asym_id']:

                    nomenclature = poly_entity['nomenclature']
                    asym_ids     = poly_entity['asym_ids']
                    auth_asym_id = poly_entity['auth_asym_id']
                    seq          = poly_entity['entity_poly_seq_one_letter_code']

                    nbr_dict[c] = BindingSiteChain(
                        sequence=seq,
                        nomenclature=nomenclature,
                        asym_ids=asym_ids,
                        auth_asym_id=auth_asym_id,
                        residues=sorted(
                            [residue for residue in nbr_residues if residue.get_parent_auth_asym_id() == c],
                            key=operator.attrgetter('seqid')
                        ))

    return BindingSite(__root__=nbr_dict)

def __getLigandResIds(ligchemid: str, struct: Structure) -> list[Residue]:
    ligandResidues: list[Residue] = list(
        filter(lambda x: x.get_resname() == ligchemid, list(struct.get_residues())))
    return ligandResidues

#※----------------------------------------------------------------------------

def struct_liglike_ids(struct_profile:RibosomeStructure) -> list[Polymer]:
    """Given an rcsb id, open the profile for the corresponding structure
    and return references to all polymers marked ligand-like"""
    pprint(struct_profile)
    polymers =  itertools.chain(struct_profile.rnas if struct_profile.rnas != None else [] , struct_profile.proteins)
    return list(filter(lambda poly : poly.ligand_like == True, polymers))

def struct_ligand_ids(pdbid: str, profile:RibosomeStructure) -> list[tuple]:
    pdbid = pdbid.upper()
    if not profile.ligands: return []
    _ = [* map(lambda x: (x.chemicalId, x.chemicalName),profile.ligands)] 
    return [ ] if len(_) < 1 else [* filter(lambda k: "ion" not in k[1].lower(), _)] 

def render_liglike_polymer(rcsb_id:str, auth_asym_id:str, structure:Structure, save:bool=False)->BindingSite:
    residues: list[Residue] = __get_polymer_residues(auth_asym_id, structure)
    binding_site_polymer: BindingSite = __get_poly_nbrs(residues, structure, auth_asym_id)

    if save:
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'POLYMER_{auth_asym_id}.json')
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)
        else:
            binding_site_polymer.json()
            with open(outfile_json, 'w') as outfile:
                json.dump(binding_site_polymer, outfile, indent=4)

    return binding_site_polymer

def render_ligand(rcsb_id:str,chemicalId:str, structure:Structure, save:bool=False)->BindingSite:
    chemicalId = chemicalId.upper()
    residues: list[Residue] = __getLigandResIds(chemicalId, structure)
    binding_site_ligand: BindingSite   = __get_ligand_nbrs(residues, structure)

    if save:
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'LIGAND_{chemicalId}.json')
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)
        else:
            binding_site_ligand.json()
            with open(outfile_json, 'w') as outfile:
                json.dump(binding_site_ligand, outfile, indent=4)

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
