from io import StringIO
from pprint import pprint
import tempfile
import os
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
import subprocess
import typing
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from functools import reduce
from io import StringIO
import os
import subprocess
from typing import Callable, Iterator, Literal
from ribctl import MUSCLE_BIN, NCBI_TAXA_SQLITE, TAXID_ARCHAEA, TAXID_BACTERIA, TAXID_EUKARYOTA
from ribctl.lib.schema.types_ribosome import PhylogenyRank
from ete3 import NCBITaxa
import os


TAXID_BACTERIA  = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA   = 2157

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
        return ncbi.get_taxid_translator([taxid])

    @staticmethod
    def get_lineage(taxid):
        return ncbi.get_lineage(taxid)

    @staticmethod
    def rank(taxid: int) -> str:
        """Given a @taxid, return the rank of the taxid"""
        lineage = ncbi.get_lineage(taxid)
        return ncbi.get_rank(lineage)[taxid]

    @staticmethod
    def ancestor_at_rank(taxid: int, target_rank: PhylogenyRank) -> int | None:
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
                new.append(Taxid.ancestor_at_rank(taxid, level))
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

class Fasta:
    records: list[SeqRecord]

    def __init__(
        self,
        path: str | None = None,
        records: list[SeqRecord] | None = None,
        # taxid_getter: Callable[[SeqRecord], int] | None = None,
    ) -> None:
        if path is not None:
            try:
                with open(path, "r") as fasta_in:
                    self.records: list[SeqRecord] = [*SeqIO.parse(fasta_in, "fasta")]
            except FileNotFoundError:
                print(f"File not found: {path}")
            except Exception as e:
                print(f"An error occurred: {str(e)}")
                exit(-1)

        elif records is not None:
            self.records = records

        # if taxid_getter is not None:
        #     self.taxid_getter = taxid_getter

    def _yield_subset(self, predicate: Callable[[SeqRecord], bool]) -> list[SeqRecord]:
        return [*filter(predicate, self.records)]

    def pick_descendants_of_taxid(self, taxid: int) -> list[SeqRecord]:
        """Given a taxid, return the subset of records that are descendants of that taxid"""
        return self._yield_subset(
            lambda record: Taxid.is_descendant_of(taxid, int(record.id))
        )

    @staticmethod
    def write_fasta(seqrecords: list[SeqRecord], outfile: str):
        with open(outfile, "w") as fasta_out:
            SeqIO.write(seqrecords, fasta_out, "fasta")

    @staticmethod
    def fasta_display_species(taxids: list[int]):
        taxids = Taxid.coerce_all_to_rank(taxids, "species")
        tree   = ncbi.get_topology(taxids)
        for node in tree.traverse():
            taxid = int(node.name)
            scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
            node.name = scientific_name
        print(tree.get_ascii(attributes=["name", "sci_name"]))

    def pick_taxids(self, taxids: list[str]) -> list[SeqRecord]:
        for taxid in set(taxids):
            if taxid not in self.all_taxids():
                raise Exception(
                    f"Taxid {taxid} not found in records. Violated assumption. Did you alter the fast archives recently?"
                )

        # ? ------ Filter duplicates. I discovered some duplicate sequences(and taxids) which breaks the HMM pipeline later on.
        # ? A cleaner solution would be to remove the duplicates from the fasta archives. Later.
        # TODO
        def filter_duplicates(acc, record):
            if record.id not in acc["seen"]:
                acc["seen"].add(record.id)
                acc["result"].append(record)
            return acc

        initial_accumulator = {"seen": set(), "result": []}
        filtered_records = reduce(
            filter_duplicates,
            filter(lambda record: record.id in taxids, self.records),
            initial_accumulator,
        )
        # ? -------------------------------------------------------------------------------

        return list(filtered_records["result"])

    def all_taxids(self) -> list[int]:
        """Given a fasta file, return all the taxids present in it
        With the assumption that the tax id is the id of each seq record."""
        taxids = []
        for record in self.records:
            # taxids = [*taxids, self.taxid_getter(record)]
            taxids = [*taxids, record.id]
        return taxids


def generate_consensus(records):
    # Convert the list of sequence records to a MultipleSeqAlignment object
    alignment = MultipleSeqAlignment(records)
    summary_align = AlignInfo.SummaryInfo(alignment)
    summary_align.dumb_consensus()
    return summary_align

def util__backwards_match(aligned_target: str, resid: int):
    """Returns the target-sequence index of a residue in the (aligned) target sequence
    Basically, "count back ignoring gaps" until you arrive at @resid
    """
    if resid > len(aligned_target):
        exit(
            IndexError(
                f"Passed residue with invalid index ({resid}) to back-match to target.Seqlen:{len(aligned_target)}"
            )
        )
    counter_proper = 0
    for i, char in enumerate(aligned_target):
        if i == resid:
            return counter_proper
        if char == "-":
            continue
        else:
            counter_proper += 1

def util__forwards_match(aligned_seq: str, resid: int):
    """Returns the index of a source-sequence residue in the aligned source sequence.
    Basically, "count forward including gaps until you reach @resid"
    """

    count_proper = 0
    for alignment_indx, char in enumerate(aligned_seq):
        if count_proper == resid:
            return alignment_indx
        if char == "-":
            continue
        else:
            count_proper += 1

def barr2str(bArr):
    return "".join([x.decode("utf-8") for x in bArr])

def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq

    _seq = _seq.replace("\n", "")
    seq_record = SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, "fasta")

def phylogenetic_neighborhood(
    taxids_base: list[str], taxid_target: str, n_neighbors: int = 10
) -> list[str]:
    """Given a set of taxids and a target taxid, return a list of the [n_neighbors] phylogenetically closest to the target."""

    tree = ncbi_.get_topology(list(set([*taxids_base, str(taxid_target)])))
    target_node = tree.search_nodes(name=str(taxid_target))[0]
    phylo_all_nodes = [(node.name, tree.get_distance(target_node, node)) for node in tree.traverse()]
    phylo_extant_nodes = filter(lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)
    phylo_sorted_nodes = sorted(phylo_extant_nodes, key=lambda x: x[1])

    nbr_taxids = list(map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

    if len(nbr_taxids) < n_neighbors:
        return nbr_taxids[1:]
    else:
        return nbr_taxids[1 : n_neighbors + 1]

def muscle_align_N_seq(
    seq_records: list[SeqRecord], vvv: bool = False
) -> Iterator[SeqRecord]:
    """Given a MSA of a protein class, and a fasta string of a chain, return a new MSA with the chain added to the class MSA."""

    with tempfile.NamedTemporaryFile(delete=False, mode="w") as temp_file:
        temp_filename = temp_file.name

        if vvv:
            pprint("Aligning sequences with muscle...")
            pprint([*seq_records])

        SeqIO.write([*seq_records], temp_filename, "fasta")
        muscle_cmd = [MUSCLE_BIN, "-in", temp_filename, "-quiet"]
        try:
            process = subprocess.run(muscle_cmd, stdout=subprocess.PIPE, text=True)
            if process.returncode == 0:
                muscle_out = process.stdout
                muscle_output_handle = StringIO(muscle_out)
                seq_records_a = SeqIO.parse(muscle_output_handle, "fasta")

                if vvv:
                    pprint("Done.")
                    pprint(seq_records_a)

                temp_file.close()
                os.remove(temp_filename)
                return seq_records_a
            else:
                print(
                    "{} failed with code {}".format(" ".join(muscle_cmd)),
                    process.returncode,
                )
                raise Exception("`{}` did not succeed.".format(" ".join(muscle_cmd)))
        except:
            temp_file.close()
            os.remove(temp_filename)
            raise Exception("Error running muscle.")

