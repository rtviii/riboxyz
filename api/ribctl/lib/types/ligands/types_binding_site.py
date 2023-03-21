import typing
from pydantic import BaseModel
from pydantic import parse_obj_as

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
NUCLEOTIDES = typing.Literal['A', 'T', 'C', 'G', 'U']



@dataclass(unsafe_hash=True, order=True)
class BindingSiteChain: 

      sequence     : str                = field(hash=True, compare=False)
      nomenclature : list[str]          = field(hash=True, compare=False)
      asym_ids     : list[str]          = field(hash=True, compare=False)
      auth_asym_id : str                = field(hash=True, compare=False)
      residues    : list[ ResidueLite ] = field(hash=True, compare=False)

class BindingSite:

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
            print(f"Saved  \033[91m{pathtofile}\033[0m.")

    def to_csv(self, pathtofile: str) -> None:
        k = [
            "chainname",
            "nomenclature",
            "residue_id",
            "residue_name"
        ]




        serialized = dict.fromkeys(k, [])

#TODO: Replace this with an actual polymer class.
@dataclass(unsafe_hash=True, order=True)
class __PolymerRef: 
      parent_rcsb_id          : str = field(hash=True, compare=False)
      auth_asym_id            : str = field(hash=True, compare=True)
      rcsb_pdbx_description   : str = field(hash=True, compare=False)
      entity_poly_seq_length  : int = field(hash=True, compare=False)
      entity_poly_polymer_type: str = field(hash=True, compare=False)