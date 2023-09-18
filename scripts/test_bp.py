import os
from pprint import pprint
from typing import Iterator
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio import SeqIO
import re

from ete3 import NCBITaxa
from ribctl import ASSETS
from ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.ribosome_assets import RibosomeAssets


class Fasta:
    records = list[SeqRecord]

    def __init__(self, path:str) -> None:
        try:
            with open(path, "r") as fasta_in:
                self.records = [*SeqIO.parse(fasta_in, "fasta")]
        except FileNotFoundError:
            print(f"File not found: {path}")
        except Exception as e:
            print(f"An error occurred: {str(e)}")
    def pick_taxids(self, taxids: list[str]) -> list[SeqRecord]:
        return list(filter(lambda record: record.id in taxids, self.records))
    def all_taxids(self) ->list[int]:
        taxids = []
        for record in self.records:
            taxids =[*taxids, record.id]
        return taxids


def process_fasta(input_file, output_file):
    try:
        with open(input_file, "r") as fasta_in, open(output_file, "w") as fasta_out:
            for seq_record in SeqIO.parse(fasta_in, "fasta"):
                # Extract the taxonomic ID from the description using regular expressions
                match = re.search(r'\|(\d+)\s*$', seq_record.description)
                if match:
                    taxonomic_id = match.group(1)
                else:
                    taxonomic_id = seq_record.id  # Use the original ID if no taxonomic ID is found

                # Remove dashes from the sequence
                unaligned_sequence = seq_record.seq.replace('-', '')

                # Create a new SeqRecord with the taxonomic ID as the ID and the unaligned sequence
                new_seq_record = SeqIO.SeqRecord(unaligned_sequence, id=taxonomic_id, description='')

                print("Writing new record", new_seq_record)
                # Write the modified record to the output file
                SeqIO.write(new_seq_record, fasta_out, "fasta")
    except FileNotFoundError:
        print(f"File not found: {input_file}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")

# for prot in list_ProteinClass:
#     fasta_file_path_in = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{prot}.fasta")
#     fasta_file_path_out = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{prot}_.fasta")
#     process_fasta(fasta_file_path_in, fasta_file_path_out)


def phylogenetic_neighborhood(taxids_base: list[str], taxid_target: str, n_neighbors: int = 10) -> list[str]:
    """Given a set of taxids and a target taxid, return a list of the [n_neighbors] phylogenetically closest to the target."""
    tree = NCBITaxa().get_topology(list(set([*taxids_base, str(taxid_target)])))
    target_node = tree.search_nodes(name=str(taxid_target))[0]
    phylo_all_nodes = [
        (node.name, tree.get_distance(target_node, node))
        for node in tree.traverse()]

    phylo_extant_nodes = filter(
        lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)

    phylo_sorted_nodes = sorted(phylo_extant_nodes, key=lambda x: x[1])

    nbr_taxids = list(
        map(lambda tax_phydist: tax_phydist[0], phylo_sorted_nodes))

    # the first element is the target node
    if len(nbr_taxids) < n_neighbors:
        return nbr_taxids[1:]
    else:
        return nbr_taxids[1:n_neighbors+1]

#!% ------------------------------------------------------------------------------------------------------------------------------------------------

def fasta_get_taxids(infile:str)->list[int]:
    """Given a fasta file, return all the taxids present in it
    With the assumption that the tax id is the id of each seq record."""
    taxids = []
    try:
        with open(infile, "r") as fasta_in:
            for seq_record in SeqIO.parse(fasta_in, "fasta"):
                taxids =[*taxids, seq_record.id]
                # Extract the taxonomic ID from the description using regular expressions
    except FileNotFoundError:
        print(f"File not found: {infile}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
    return taxids


rib      = RibosomeAssets('3J7Z').profile()
organism = rib.src_organism_ids[0]
prots    = RibosomeAssets('3J7Z').profile().proteins

i = 0
for protclass in list_ProteinClass:
    fasta_path = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{protclass}.fasta")
    ids = fasta_get_taxids(fasta_path)
    phylo_nbhd=  phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism), n_neighbors=10)
    Fasta(fasta_path).pick_taxids(phylo_nbhd)

    i+=1
    if i > 5:
        exit()



fasta_path = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"bL9.fasta")
record = [*SeqIO.parse(fasta_path, "fasta")]




# for prot in prots:


    









# extract sequences 
# store (create an asset class in __init__), possibly compress
# create method to pick phyl. nbhd, create an MSA, create an HMM
# ->evaluate as before