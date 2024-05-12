from typing_extensions import Literal
from Bio.Align import MultipleSeqAlignment, AlignInfo
import typing
from Bio.SeqRecord import SeqRecord
import subprocess
from pydantic import BaseModel
from ribctl import MUSCLE_BIN, NCBI_TAXA_SQLITE
import os
from ete3 import NCBITaxa

# ? ---------------- [ Phylogeny] ----------------

TAXID_BACTERIA  = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA   = 2157
PhylogenyRank = Literal["superkingdom", "phylum", "class", "order", "family", "genus", "species", "strain"]
ncbi          = NCBITaxa(dbfile=NCBI_TAXA_SQLITE)

class Taxid:
    @staticmethod
    def is_descendant_of(parent_taxid: int, target_taxid: int) -> bool:
        lineage = ncbi.get_lineage(target_taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        return False if parent_taxid not in lineage else True

    @staticmethod
    def get_name(taxid):
        return list(ncbi.get_taxid_translator([taxid]).values())[0]

    @staticmethod
    def get_lineage(taxid)->list[int]:
        """Return ncbi lineage, except filter out the ranks that are not among the @PhylogenyRank."""
        return list(filter(lambda x: Taxid.rank(x) in typing.get_args(PhylogenyRank), ncbi.get_lineage(taxid) ) )

    @staticmethod
    def rank(taxid: int) -> PhylogenyRank:
        """Given a @taxid, return the rank of the taxid"""
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage)[taxid]

    @staticmethod
    def coerce_to_rank(taxid: int, target_rank: PhylogenyRank) -> int | None:
        """Given a @taxid and a @rank, return the taxid of the first ancestor of @taxid that is at @rank"""
        lineage = ncbi.get_lineage(taxid)
        if lineage is None:
            raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
        for item in lineage:
            rank = ncbi.get_rank([item])[item]
            if rank == target_rank:
                return item

        raise IndexError("Taxid {} does not have a {} level".format(taxid, target_rank))

    @staticmethod
    def superkingdom(
        taxid: int,
    ) -> typing.Literal["bacteria", "eukaryota", "archaea", "virus"]:
        match (
            Taxid.is_descendant_of(TAXID_EUKARYOTA, taxid),
            Taxid.is_descendant_of(TAXID_BACTERIA, taxid),
            Taxid.is_descendant_of(TAXID_ARCHAEA, taxid),
        ):
            case (False, False, True):
                return "archaea"
            case (False, True, False):
                return "bacteria"
            case (True, False, False):
                return "eukaryota"
            case (False, False, False):
                print("Probably a virus")
                return "virus"
            case _:
                raise ValueError(
                    "Taxid {} is not a descendant of any of the three domains".format(
                        taxid
                    )
                )

    @staticmethod
    def coerce_all_to_rank(taxids: list[int], level: PhylogenyRank) -> list[int]:
        """Given a list of taxids, return a list of the same taxids but coerced to the given lineage level(rank)."""
        new = []
        for taxid in taxids:
            try:
                new.append(Taxid.coerce_to_rank(taxid, level))
            except Exception as e:
                print(e)
                raise Exception(e)
        return new

    @staticmethod
    def get_descendants_of(parent: int, targets: list[int]):
        """Given a @parent taxid and a list of @taxids, return the subset of @taxids that are descendants of @parent"""
        descendants = set()
        for tax_id in targets:
            lineage = ncbi.get_lineage(tax_id)
            if lineage == None:
                raise LookupError("Lineage is None. Check if taxid is NCBI-valid.")
            if parent in lineage:
                descendants.add(tax_id)
        return descendants

class PhylogenyNode(BaseModel):

    @staticmethod
    def from_taxid(taxid:int):
        return PhylogenyNode(ncbi_tax_id=taxid, scientific_name=Taxid.get_name(taxid), rank=Taxid.rank(taxid))

    def __hash__(self) -> int:
        return self.ncbi_tax_id

    def get_lineage(self) -> list[int]:
        return Taxid.get_lineage(self.ncbi_tax_id)

    def get_lineage_nodes(self) -> list[int]:
        _ = []
        for taxid in self.get_lineage():
            if Taxid.rank(taxid) not in typing.get_args(PhylogenyRank):
               continue
            _.append(self.from_taxid(taxid))
        return _

    ncbi_tax_id    : int
    scientific_name: str
    rank           : PhylogenyRank


# ? ----------------------------------------------{ Subcomponent Types }------------------------------------------------