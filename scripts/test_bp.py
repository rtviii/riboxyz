import os
from pprint import pprint
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio import SeqIO
import re

from ribctl import ASSETS
from ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.ribosome_assets import RibosomeAssets


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



rib = RibosomeAssets('3J7Z').profile()
organism = rib.src_organism_ids[0]
prots = RibosomeAssets('3J7Z').profile().proteins

i = 0
for protclass in list_ProteinClass:
    ids = fasta_get_taxids(os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{protclass}.fasta"))
    
    i+=1
    if i > 5:
        exit()
    # print(protclass)



# for prot in prots:


    









# extract sequences 
# store (create an asset class in __init__), possibly compress
# create method to pick phyl. nbhd, create an MSA, create an HMM
# ->evaluate as before