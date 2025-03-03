from io import StringIO
import sys
sys.path.append('/home/rtviii/dev/riboxyz')
from pprint import pprint
import tempfile
import os
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment, AlignInfo
import subprocess
import typing
from Bio.SeqRecord import SeqRecord
from functools import reduce
import os
import subprocess
from typing import Callable, Iterator
from ribctl import ASSETS, MUSCLE_BIN
import os
from ribctl.lib.libtax import  Taxid, get_ncbi
from ribctl.lib.types.polymer import CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, MitochondrialProteinClass, PolymerClass, PolynucleotideClass
import os
from typing import Dict, List, Set, Tuple
from Bio import SeqIO
from collections import defaultdict
from ete3 import NCBITaxa, Tree
import random
from rich.console import Console
from rich.tree import Tree as RichTree
from rich import print as rprint
from rich.panel import Panel
from rich.table import Table

from ribctl import ASSETS, NCBI_TAXA_SQLITE
from ribctl.lib.libtax import Taxid


class Fasta :
    records: list[SeqRecord]
    def __init__( self, path: str | None = None, records: list[SeqRecord] | None = None ) -> None:
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

    @staticmethod
    def poly_class_all_seq(candidate_class:PolymerClass):
        
        if candidate_class in CytosolicProteinClass:
            fasta_path = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")
        elif candidate_class in MitochondrialProteinClass:
            fasta_path = os.path.join(ASSETS["fasta_proteins_mitochondrial"], f"{candidate_class.value}.fasta")
        elif candidate_class in PolynucleotideClass:
            fasta_path = os.path.join(ASSETS["fasta_rna"], f"{candidate_class.value}.fasta")
        elif candidate_class in ElongationFactorClass:
            fasta_path = os.path.join(ASSETS["fasta_factors_elongation"], f"{candidate_class.value}.fasta")
        elif candidate_class in InitiationFactorClass:
            fasta_path = os.path.join(ASSETS["fasta_factors_initiation"], f"{candidate_class.value}.fasta")
        else:
            raise KeyError(f"Class {candidate_class} not found in any of the fasta archives. Something went terribly wrong.")

        fasta_path = os.path.join(ASSETS["fasta_proteins_cytosolic"], f"{candidate_class.value}.fasta")
        return Fasta(fasta_path)

    def _yield_subset(self, predicate: Callable[[SeqRecord], bool]) -> list[SeqRecord]:
        return [*filter(predicate, self.records)]

    def pick_descendants_of_taxid(self, taxid: int) -> list[SeqRecord]:
        """Given a taxid, return the subset of records that are descendants of that taxid"""
        return self._yield_subset( lambda record: Taxid.is_descendant_of(taxid, int(record.id)))

    @staticmethod
    def write_fasta(seqrecords: list[SeqRecord], outfile: str):
        with open(outfile, "w") as fasta_out:
            SeqIO.write(seqrecords, fasta_out, "fasta")

    @staticmethod
    def fasta_display_species(taxids: list[int]):
        taxids = Taxid.coerce_all_to_rank(taxids, "species")
        ncbi = get_ncbi()
        tree   = ncbi.get_topology(taxids)

        for node in tree.traverse():
            taxid           = int(node.name)
            scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
            node.name       = scientific_name

    def pick_taxids(self, _taxids_int: list[int]) -> list[SeqRecord]:
         
        taxids = list(map(lambda x: str(x), _taxids_int))
        for taxid in set(taxids):
            if taxid not in self.all_taxids(return_as='str'):
                raise Exception( f"Taxid {taxid} not found in records. Violated assumption. Did you alter the fast archives recently? There might be a cached HMM Scanner with taxid that is no longer present. " )

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

    def all_taxids(self, return_as:typing.Literal["int", "str"]="str") -> list[int] | list[str]:
        """Given a fasta file, return all the taxids present in it
        With the assumption that the tax id is the id of each seq record."""
        taxids = []
        for record in self.records:
            # taxids = [*taxids, self.taxid_getter(record)]
            taxids = [*taxids, record.id]
        if return_as == "str":
            return taxids
        elif return_as == "int":
            return list(map(lambda _: int(_), taxids))
        else:
            raise Exception("Invalid type passed to all_taxids")
    
def generate_consensus(records):
    alignment     = MultipleSeqAlignment(records)
    summary_align = AlignInfo.SummaryInfo(alignment)
    summary_align.dumb_consensus()
    return summary_align

def barr2str(bArr):
    return "".join([x.decode("utf-8") for x in bArr])

def seq_to_fasta(rcsb_id: str, _seq: str, outfile: str):
    from Bio.Seq import Seq
    _seq = _seq.replace("\n", "")
    seq_record = SeqRecord(Seq(_seq).upper())
    seq_record.id = seq_record.description = rcsb_id
    SeqIO.write(seq_record, outfile, "fasta")


def muscle_align_N_seq( seq_records: list[SeqRecord], vvv: bool = False ) -> Iterator[SeqRecord]:
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
                muscle_out           = process.stdout
                muscle_output_handle = StringIO(muscle_out)
                seq_records_a        = SeqIO.parse(muscle_output_handle, "fasta")

                if vvv:
                    pprint("Done.")
                    pprint(seq_records_a)

                temp_file.close()
                os.remove(temp_filename)
                return seq_records_a
            else:
                print( "{} failed with code {}".format(" ".join(muscle_cmd)), process.returncode, )
                raise Exception("`{}` did not succeed.".format(" ".join(muscle_cmd)))
        except:
            temp_file.close()
            os.remove(temp_filename)
            raise Exception("Error running muscle.")


class FastaBalancer:
    def __init__(self, fasta: Fasta):
        """Initialize with a Fasta instance containing sequences with taxids as headers."""
        self.fasta = fasta
        self.ncbi = NCBITaxa(dbfile=NCBI_TAXA_SQLITE)
        self.sequences = self._load_sequences()
        self.console = Console()

    def _load_sequences(self) -> Dict[int, str]:
        """Load sequences from Fasta instance."""
        return {int(record.id): str(record.seq) for record in self.fasta.records}

    def get_taxonomic_hierarchy(self) -> Dict[str, Dict[str, List[int]]]:
        """Build complete taxonomic hierarchy of sequences."""
        hierarchy = defaultdict(lambda: defaultdict(list))

        for taxid in self.sequences.keys():
            try:
                kingdom = Taxid.superkingdom(taxid)
                lineage = Taxid.get_lineage(taxid)
                phylum_id = next(
                    (tid for tid in lineage if Taxid.rank(tid) == "phylum"), None
                )

                if phylum_id:
                    phylum_name = Taxid.get_name(phylum_id)
                    hierarchy[kingdom][phylum_name].append(taxid)
            except Exception as e:
                print(f"Warning: Could not process taxid {taxid}: {str(e)}")

        return hierarchy

    def balance_dataset_sparse(self, target_size: int) -> Fasta:
        """
        Balance dataset by maximizing phylogenetic diversity for a given target size.
        Returns a new Fasta instance with only selected sequences.

        Args:
            target_size: Desired total number of sequences
        """
        original_taxids = set(self.sequences.keys())
        hierarchy = self.get_taxonomic_hierarchy()

        self.console.print(
            Panel(
                f"[bold]Starting sparse dataset balancing[/bold]\nTarget total sequences: {target_size}"
            )
        )

        selected_taxids = set()

        total_kingdoms = len(hierarchy)
        total_phyla = sum(len(phyla) for phyla in hierarchy.values())

        if target_size <= total_kingdoms:
            for kingdom in list(hierarchy.keys())[:target_size]:
                phylum_sizes = [
                    (p, len(taxa)) for p, taxa in hierarchy[kingdom].items()
                ]
                largest_phylum = max(phylum_sizes, key=lambda x: x[1])[0]
                selected = random.sample(hierarchy[kingdom][largest_phylum], 1)
                selected_taxids.update(selected)
        else:
            seqs_per_kingdom = target_size // total_kingdoms
            remaining = target_size % total_kingdoms

            for kingdom, phyla in hierarchy.items():
                kingdom_target = seqs_per_kingdom + (1 if remaining > 0 else 0)
                remaining = max(0, remaining - 1)

                if kingdom_target == 0:
                    continue

                n_phyla = len(phyla)
                seqs_per_phylum = max(1, kingdom_target // n_phyla)
                phylum_remainder = kingdom_target % n_phyla

                self.console.print(f"\n[bold]{kingdom}[/bold]:")
                self.console.print(
                    f"Found {n_phyla} phyla, allocating {kingdom_target} sequences"
                )

                for phylum_name, phylum_taxids in phyla.items():
                    phylum_target = min(
                        seqs_per_phylum + (1 if phylum_remainder > 0 else 0),
                        len(phylum_taxids),
                    )
                    phylum_remainder = max(0, phylum_remainder - 1)

                    selected = random.sample(phylum_taxids, phylum_target)
                    selected_taxids.update(selected)
                    self.console.print(
                        f"  • {phylum_name}: selected {len(selected)}/{len(phylum_taxids)} sequences"
                    )

        self.print_filtering_summary(original_taxids, selected_taxids)

        # Create new Fasta instance with only selected sequences

        selected_records = [record for record in self.fasta.records 
                           if int(record.id) in selected_taxids]
        return Fasta(records=selected_records)

    def print_taxonomic_tree(self, taxids: List[int], title: str):
        """Print a hierarchical view of the taxonomic distribution."""
        rich_tree = RichTree(f"[bold blue]{title}[/bold blue]")

        kingdom_groups = defaultdict(list)
        for taxid in taxids:
            try:
                kingdom = Taxid.superkingdom(taxid)
                kingdom_groups[kingdom].append(taxid)
            except Exception:
                continue

        for kingdom, kingdom_taxids in kingdom_groups.items():
            kingdom_branch = rich_tree.add(
                f"[bold green]{kingdom}[/bold green] ({len(kingdom_taxids)} sequences)"
            )

            phylum_groups = defaultdict(list)
            for taxid in kingdom_taxids:
                try:
                    lineage = Taxid.get_lineage(taxid)
                    phylum_id = next(
                        (tid for tid in lineage if Taxid.rank(tid) == "phylum"), None
                    )
                    if phylum_id:
                        phylum_groups[Taxid.get_name(phylum_id)].append(taxid)
                except Exception:
                    continue

            for phylum, phylum_taxids in phylum_groups.items():
                phylum_branch = kingdom_branch.add(
                    f"[yellow]{phylum}[/yellow] ({len(phylum_taxids)} sequences)"
                )

        self.console.print(rich_tree)

    def print_filtering_summary(
        self, original_taxids: Set[int], selected_taxids: Set[int]
    ):
        """Print detailed summary of filtering process."""
        filtered_taxids = original_taxids - selected_taxids

        table = Table(title="Filtering Summary")
        table.add_column("Category", style="cyan")
        table.add_column("Count", justify="right", style="green")
        table.add_column("Percentage", justify="right", style="yellow")

        total = len(original_taxids)
        selected = len(selected_taxids)
        filtered = len(filtered_taxids)

        table.add_row("Original sequences", str(total), "100%")
        table.add_row(
            "Selected sequences", str(selected), f"{(selected/total)*100:.1f}%"
        )
        table.add_row(
            "Filtered sequences", str(filtered), f"{(filtered/total)*100:.1f}%"
        )

        self.console.print("\n")
        self.console.print(table)

        self.console.print("\n[bold]Taxonomic Distribution Before Filtering:[/bold]")
        self.print_taxonomic_tree(list(original_taxids), "Original Dataset")

        self.console.print("\n[bold]Taxonomic Distribution After Filtering:[/bold]")
        self.print_taxonomic_tree(list(selected_taxids), "Balanced Dataset")

# if __name__ == "__main__":
#     fasta          = Fasta(os.path.join(ASSETS["fasta_proteins_cytosolic"], "uL4.fasta"))
#     balancer       = FastaBalancer(fasta)
#     filtered_fasta = balancer.balance_dataset_sparse(target_size=10)
