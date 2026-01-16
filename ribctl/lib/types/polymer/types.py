from typing import Union
from ribctl.lib.types.polymer.base import (
    CytosolicProteinClass,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    MitochondrialProteinClass,
    MitochondrialRNAClass,
    tRNA,
)

ProteinClass         = Union[MitochondrialProteinClass, CytosolicProteinClass]
LifecycleFactorClass = Union[ElongationFactorClass, InitiationFactorClass]
PolypeptideClass     = Union[LifecycleFactorClass, ProteinClass]
PolynucleotideClass  = Union[CytosolicRNAClass, MitochondrialRNAClass, tRNA]
PolymerClass         = Union[PolynucleotideClass, PolypeptideClass]
