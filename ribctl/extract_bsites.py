import dataclasses
import json
import os,sys
from pprint import pprint
sys.path.append(os.environ.get("PYMOL_PATH")) 
RIBETL_DATA = str(os.environ.get('RIBETL_DATA'))
from typing import Dict, List
import operator
from Bio.PDB import MMCIF2Dict
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Residue import Residue
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
import argparse
import itertools
from dataclasses import dataclass, field
import itertools
from utils import  open_structure, struct_path
flatten = itertools.chain.from_iterable

def __get_dict(path:str,)->dict:
	return MMCIF2Dict.MMCIF2Dict(path)

AMINO_ACIDS = {
    "ALA": 0,
    'ARG': 1,
    'ASN': 0,
    'ASP': -1,
    'CYS': 0,
    'GLN': 0,
    'GLU': -1,
    'GLY': 0,
    'HIS': 0,
    'ILE': 0,
    'LEU': 0,
    'LYS': 1,
    'MET': 0,
    'PHE': 0,
    'PRO': 0,
    'SER': 0,
    'THR': 0,
    'TRP': 0,
    'TYR': 0,
    'VAL': 0,
    'SEC': 0,
    'PYL': 0
    }
NUCLEOTIDES = ['A', 'T', 'C', 'G', 'U']


def __open_ligand(pdbid: str, ligid: str, ligpath: str |None = None):
    pdbid = pdbid.upper()
    ligid = ligid.upper()

    if ligpath == None:
        ligpath = os.path.join(RIBETL_DATA, pdbid, f'LIGAND_{ligid}.json')
    with open(ligpath, 'rb') as infile:
        data = json.load(infile)
    return data

@dataclass(unsafe_hash=True, order=True)
class ResidueLite:

    residue_name        : str = field(hash=True, compare=False)
    residue_id          : int = field(hash=True, compare=True)
    parent_auth_asym_id : str  = field(hash=True, compare=False)

    @staticmethod
    def res2reslite(r: Residue):
        biopy_id_tuple = r.get_full_id()
        parent_chain   = biopy_id_tuple[2]
        resname        = r.resname
        resid          = r.id[1]

        return ResidueLite(resname, resid, parent_chain)

@dataclass(unsafe_hash=True, order=True)
class BindingSiteChain: 
      sequence        : str = field(hash=True, compare=False)
      nomenclature    : list[str]= field(hash=True, compare=False)
      asym_ids        : list[str]= field(hash=True, compare=False)
      auth_asym_id    : str= field(hash=True, compare=False)
      residues: list[ ResidueLite ]         = field(hash=True, compare=False)

class __BindingSite:

    def __init__(self, data: Dict[str, BindingSiteChain]) -> None:
        self.data: Dict[str, BindingSiteChain] = data

    def __getitem__(self, chainkey:str)->BindingSiteChain:
        if chainkey not in self.data:
            raise KeyError(f"Chain {chainkey} not found in the binding site.")
        return self.data[chainkey]

    def to_json(self, pathtofile: str) -> None:
        with open(pathtofile, 'w') as outf:
            serialized = {}
            for x in self.data.items():
                serialized.update({x[0]: dataclasses.asdict(x[1])})
            json.dump(serialized, outf)
            print(f"\033[91mSaved  {pathtofile} successfuly.\033[0m")

    def to_csv(self, pathtofile: str) -> None:
        k = [
            "chainname",
            "nomenclature",
            "residue_id",
            "residue_name"
        ]
        serialized = dict.fromkeys(k, [])

@dataclass(unsafe_hash=True, order=True)
class __PolymerRef: 
      parent_rcsb_id          : str = field(hash=True, compare=False)
      auth_asym_id            : str = field(hash=True, compare=True)
      rcsb_pdbx_description   : str = field(hash=True, compare=False)
      entity_poly_seq_length  : int = field(hash=True, compare=False)
      entity_poly_polymer_type: str = field(hash=True, compare=False)

def __get_polymer_residues(auth_asym_id: str, struct: Structure) -> List[Residue]:
    c: Chain = struct[0][auth_asym_id]
    return [*c.get_residues()]

def __get_poly_nbrs(
      residues      : List[Residue],
      struct        : Structure,
      auth_asym_id  : str
    ) -> __BindingSite: 

    """KDTree search the neighbors of a given list of residues(which constitue a ligand)
    and return unique having tagged them with a ban identifier proteins within 5 angstrom of these residues. """

    pdbid         = struct.get_id().upper()
    ns            = NeighborSearch(list(struct.get_atoms()))
    nbr_residues  = []
    parent_strand = residues[0].get_parent().id if len(residues) > 0 else ...

    for poly_res in residues:
        for atom in poly_res.child_list:
            nbr_residues.extend(ns.search(atom.get_coord(), 10, level='R'))

    nbr_residues = list(set([* map(ResidueLite.res2reslite, nbr_residues)]))
    nbr_residues = list(filter(lambda res: res.parent_auth_asym_id != parent_strand and res.residue_name in [*AMINO_ACIDS.keys(), *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _: _.parent_auth_asym_id, nbr_residues)))

    with open(struct_path(pdbid, 'json'), 'rb') as strfile:
        profile       = json.load(strfile)
        poly_entities = [*profile['proteins'], *profile['rnas']]

    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity['auth_asym_id']:
                    nomenclature = poly_entity['nomenclature']
                    asym_id      = poly_entity['asym_ids']
                    auth_asym_id = poly_entity['auth_asym_id']
                    seq          = poly_entity['entity_poly_seq_one_letter_code']

        nbr_dict[c] = BindingSiteChain(
            seq,
            nomenclature,
            asym_id,
            auth_asym_id,
            sorted([residue for residue in nbr_residues if residue.parent_auth_asym_id == c], key=operator.attrgetter('residue_id')))


    return __BindingSite(nbr_dict)

def __get_ligand_nbrs(
      ligand_residues: List[Residue],
      struct         : Structure,
      chemicalId: str
    ) -> __BindingSite : 
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
    nbr_residues = list(set([* map(ResidueLite.res2reslite, nbr_residues)]))

    # Filter the ligand itself, water and other special residues
    nbr_residues = list(filter(lambda resl: resl.residue_name in [*AMINO_ACIDS.keys(),  *NUCLEOTIDES], nbr_residues))
    nbr_dict     = {}
    chain_names  = list(set(map(lambda _:  _.parent_auth_asym_id, nbr_residues)))



    for c in chain_names:
        for poly_entity in poly_entities:
            if c == poly_entity['auth_asym_id']:

                    nomenclature = poly_entity['nomenclature']
                    asym_ids     = poly_entity['asym_ids']
                    auth_asym_id = poly_entity['auth_asym_id']
                    seq          = poly_entity['entity_poly_seq_one_letter_code']

        nbr_dict[c] = BindingSiteChain(
            seq,
            nomenclature,
            asym_ids,
            auth_asym_id,
            sorted([residue for residue in nbr_residues if residue.parent_auth_asym_id == c], key=operator.attrgetter('residue_id')))


    return __BindingSite(nbr_dict)

def __getLigandResIds(ligchemid: str, struct: Structure) -> List[Residue]:
    ligandResidues: List[Residue] = list(
        filter(lambda x: x.get_resname() == ligchemid, list(struct.get_residues())))
    return ligandResidues

#※----------------------------------------------------------------------------

def get_liglike_polymers(struct_profile:dict) -> List[__PolymerRef]:
    """Given an rcsb id, open the profile for the corresponding structure
    and return references to all polymers marked ligand-like"""
    liglike = []
    for i in [*struct_profile['rnas'], *struct_profile['proteins']]:
        if 'ligand_like' not in i.keys(): raise KeyError("ligand_like key not found in struct json_profile. Perhaps the semantics have changed.")
        if i['ligand_like'] == True:
            liglike =  [*liglike,
                       __PolymerRef(
                           i['parent_rcsb_id'],
                           i['auth_asym_id'],
                           i['rcsb_pdbx_description'],
                           i['entity_poly_seq_length'],
                           i['entity_poly_polymer_type'],
                       )]
    return liglike

def get_ligands(pdbid: str, profile:dict) -> List[tuple]:
    pdbid = pdbid.upper()
    if not profile['ligands']: return []

    _ = [* map(lambda x: (x['chemicalId'], x['chemicalName']),profile['ligands'])] 
    return [ ] if len(_) < 1 else [* filter(lambda k: "ion" not in k[1].lower(), _)] 

def render_liglike_polymer(rcsb_id:str, auth_asym_id:str, structure:Structure, save:bool=False)->__BindingSite:
    residues: list[Residue] = __get_polymer_residues(auth_asym_id, structure)
    binding_site_polymer: __BindingSite = __get_poly_nbrs(residues, structure, auth_asym_id)

    if args.save:
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'POLYMER_{auth_asym_id}.json')
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)
        else:
            binding_site_polymer.to_json(outfile_json)

    return binding_site_polymer

def render_ligand(rcsb_id:str,chemicalId:str, structure:Structure, save:bool=False)->__BindingSite:
    chemicalId = chemicalId.upper()
    residues: list[Residue] = __getLigandResIds(chemicalId, structure)
    binding_site_ligand: __BindingSite   = __get_ligand_nbrs(residues, structure, chemicalId)

    if args.save:
        outfile_json = os.path.join(RIBETL_DATA, rcsb_id.upper(), f'LIGAND_{chemicalId}.json')
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)
        else:
            binding_site_ligand.to_json(outfile_json)
    return binding_site_ligand

if __name__ == "__main__":


    parser = argparse. ArgumentParser(description='Split structure into constituent polymers and inject new nomencalture into the .cif file')
    parser.add_argument ('-s'     , '--structure', type=   str   , required=True                                                            )
    parser.add_argument ('--save' ,action          ='store_true'                                                    )
    
    args  = parser.parse_args()
    PDBID = args.structure.upper()


    _structure_cif_handle :Structure = open_structure(PDBID,'cif')  # type: ignore
    struct_profile_handle:dict       = open_structure(PDBID,'json')  # type: ignore

    liglike_polys = get_liglike_polymers(struct_profile_handle)
    ligands       = get_ligands(PDBID, struct_profile_handle)


    pprint(ligands)
    pprint(liglike_polys)

    for polyref in liglike_polys:
        render_liglike_polymer(polyref.parent_rcsb_id, polyref.auth_asym_id, _structure_cif_handle, args.save)

    for l in ligands:
        print(f"Ligand {l[0]}")
        outfile_json = os.path.join(RIBETL_DATA, PDBID.upper(), f'LIGAND_{l[0].upper()}.json')
        if (os.path.isfile(outfile_json)):
            print("Exists already: ", outfile_json)
        residues: List[Residue] = __getLigandResIds(l[0].upper(), _structure_cif_handle)
        bsl      : __BindingSite   = __get_ligand_nbrs(residues, _structure_cif_handle, l[0].upper())
        if args.save:
            bsl.to_json(outfile_json)

