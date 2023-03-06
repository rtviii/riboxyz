from ninja import Schema
from api.ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import Protein, ProteinClass
"""This file documents the possible requests that the API can receive."""


# interface proteins_number
class ProteinsNumber(Schema):
    pass

class BanclassAnnotation(Schema):
    banclass: ProteinClass


class NomclassVisualize(Schema):
    ban: str


class BanclassMetadata(Schema):
    family: str
    subunit: str


class IndividualLigand(Schema):
    chemId: str


class StructsByProteins(Schema):
    proteins: str


class RnasByStruct(Schema):
    pass


class LigandsByStruct(Schema):
    pass


class RnaClass(Schema):
    rna_class: RNAClass


class Structure(Schema):
    pdbid: str


class Homologs(Schema):
    banName: Protein


class AllStructures(Schema):
    pass


class ListNomClasses(Schema):
    pass


class MembersOfProteinClass(Schema):
    banName: Protein


class AllLigandlike(Schema):
    pass


class AllLigands(Schema):
    pass


class LigandNbhd(Schema): 
      src_struct        : str
      ligandlike_id     : str
      is_polymer        : bool
