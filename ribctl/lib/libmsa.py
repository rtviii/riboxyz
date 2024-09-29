from io import StringIO
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
from ribctl.lib.libtax import  Taxid, ncbi
from ribctl.lib.schema.types_ribosome import CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, MitochondrialProteinClass, PolymerClass, PolynucleotideClass


class Fasta:
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
        tree   = ncbi.get_topology(taxids)

        for node in tree.traverse():
            taxid           = int(node.name)
            scientific_name = ncbi.get_taxid_translator([taxid]).get(taxid, "Unknown")
            node.name       = scientific_name
        print(tree.get_ascii(attributes=["name", "sci_name"]))

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

def phylogenetic_neighborhood( _taxids_base: list[int], taxid_target: str, n_neighbors: int = 10 ) -> list[int]:
    """Given a set of taxids and a target taxid, return a list of the [n_neighbors] phylogenetically closest to the target."""

    taxids_base        = list(map(lambda x: str(x),_taxids_base)) # Ensure all taxids are strings because that's what ete's ncbi expects

    tree               = ncbi.get_topology(list(set([*taxids_base, str(taxid_target)])))
    target_node        = tree.search_nodes(name=str(taxid_target))[0]
    phylo_all_nodes    = [(other_node.name, tree.get_distance(target_node, other_node)) for other_node in tree.traverse()]

    phylo_extant_nodes = filter(lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)
    phylo_sorted_nodes = sorted(phylo_extant_nodes, key=lambda x: x[1]) # Sort by phylogenetic distance (tuples are (taxid, phylo_dist) ex. ('4932', 3.0))
    nbr_taxids         = list(map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

    if len(nbr_taxids) < n_neighbors:
        return list(map(lambda x : int(x),  nbr_taxids[1:] ))
    else:
        return list(map(lambda x : int(x), nbr_taxids[1 : n_neighbors + 1]))

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








