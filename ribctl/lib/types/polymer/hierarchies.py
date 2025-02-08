from enum import Enum
from typing import  Union, Type, Iterator, TypeVar, Generic

from ribctl.lib.types.polymer.base import CytosolicProteinClass, CytosolicRNAClass, ElongationFactorClass, InitiationFactorClass, MitochondrialProteinClass, MitochondrialRNAClass, tRNA

E = TypeVar("E", bound=Enum)


class PolymerHierarchy(Generic[E]):
    """
    Manages a hierarchical collection of related Enum classes.
    Provides type-safe iteration, containment checking, and classification.
    """

    def __init__(self, *enum_classes: Type[E], name: None | str = None):
        self.enum_classes = enum_classes
        self.name = name or "+".join(cls.__name__ for cls in enum_classes)

        # Pre-compute all valid values for faster lookup
        self._all_values = {
            member.value: (cls, member) for cls in enum_classes for member in cls
        }

    def __iter__(self) -> Iterator[E]:
        """Allows iteration over all members across all enum classes."""
        return (member for cls in self.enum_classes for member in cls)

    def __contains__(self, item: Union[str, E]) -> bool:
        """Enables 'in' operator for both string values and enum members."""
        if isinstance(item, str):
            return item in self._all_values
        return item in self._all_values.values()

    def get_type(self, value: str) -> Type[E] | None:
        """Returns the enum class type for a given value."""

        if value in self._all_values:
            return self._all_values[value][0]
        return None

    def get_member(self, value: str) -> E | None:
        """Returns the enum member for a given value."""
        if value in self._all_values:
            return self._all_values[value][1]
        return None

    def is_of_type(self, value: Union[str, E], enum_type: Type[E]) -> bool:
        """Checks if a value belongs to a specific enum type."""
        if isinstance(value, str):
            cls = self.get_type(value)
            return cls == enum_type if cls else False
        return isinstance(value, enum_type)


Proteins         = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, name="Proteins" )
LifecycleFactors = PolymerHierarchy( ElongationFactorClass, InitiationFactorClass, name="LifecycleFactors" )
Polypeptides     = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, name="Polypeptides", )
Polynucleotides  = PolymerHierarchy( CytosolicRNAClass, MitochondrialRNAClass, tRNA, name="Polynucleotides" )
Polymers         = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, CytosolicRNAClass, MitochondrialRNAClass, tRNA, name="Polymers" )