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
from ribctl.lib.libmsa import Fasta
from ribctl.lib.libtax import Taxid


class SequenceBalancer:
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

    # [Previous visualization methods remain unchanged]
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


if __name__ == "__main__":

    fasta          = Fasta(os.path.join(ASSETS["fasta_rna"], "16SrRNA.fasta"))
    balancer       = SequenceBalancer(fasta)
    filtered_fasta = balancer.balance_dataset_sparse(target_size=10)
