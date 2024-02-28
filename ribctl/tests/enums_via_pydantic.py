from enum import Enum, IntEnum
from typing import Union
from pydantic import BaseModel, ValidationError
from ribctl.lib.enumunion import enum_union

class MitochondrialRNAClass(Enum):
    mt_rRNA_12S = "mt12SrRNA" # mitochondrial
    mt_rRNA_16S = "mt16SrRNA" # mitochondrial

class CytosolicRNAClass(Enum):

    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA" #  c-bacterial or mitochondrial



PolynucleotideClass = Union[MitochondrialRNAClass, CytosolicRNAClass]
Polynucleotide      = enum_union( MitochondrialRNAClass, CytosolicRNAClass)





members=  [*map( lambda x: x.value, [*Polynucleotide] )]
print(members)
x = "5SrRNA" in members
print(x)

# just keeep enum union for now.



