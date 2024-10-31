

from ribctl.etl.etl_assets_ops import RibosomeOps
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.Residue import Residue
from Bio.PDB.MMCIFParser import FastMMCIFParser


struct:Structure = RibosomeOps('4U3M').biopython_structure()
chain :Chain     = struct[0]['1']

[ print(_) for _ in chain.child_list]

