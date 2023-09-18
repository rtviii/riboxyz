from functools import reduce
from io import StringIO
from itertools import tee
import os
from pprint import pprint
import subprocess
from tempfile import NamedTemporaryFile
from typing import Iterator
from Bio.Align import MultipleSeqAlignment,Seq, SeqRecord
from Bio.Align.Applications import MuscleCommandline
from Bio import SeqIO, AlignIO, pairwise2
import re
from ete3 import NCBITaxa
import pyhmmer
from ribctl import ASSETS, MUSCLE_BIN
from ribctl.lib.types.types_poly_nonpoly_ligand import list_ProteinClass
from ribctl.lib.types.types_ribosome import ProteinClass
from ribctl.ribosome_assets import RibosomeAssets
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence, SequenceFile, SequenceBlock, TextSequenceBlock
from pyhmmer.plan7 import Pipeline, HMM


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
        
        for taxid in set(taxids): 
            if taxid not in self.all_taxids():
                raise Exception(f"Taxid {taxid} not found in records. Violated assumption. Did you alter the fast archives recently?")

        #? ------ Filter duplicates. I discovered some duplicate sequences(and taxids) which breaks the HMM pipeline later on.
        #? A cleaner solution would be to remove the duplicates from the fasta archives. Later.
        # TODO
        def filter_duplicates(acc, record):
            if record.id not in acc['seen']:
                acc['seen'].add(record.id)
                acc['result'].append(record)
            return acc

        initial_accumulator = {'seen': set(), 'result': []}
        filtered_records = reduce(filter_duplicates, filter(lambda record: record.id in taxids, self.records), initial_accumulator)
        #? -------------------------------------------------------------------------------

        return list(filtered_records['result'])

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
    phylo_all_nodes = [(node.name, tree.get_distance(target_node, node)) for node in tree.traverse()]

    phylo_extant_nodes = filter(lambda taxid: taxid[0] in taxids_base, phylo_all_nodes)

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
                taxids = [*taxids, seq_record.id]
                # Extract the taxonomic ID from the description using regular expressions
    except FileNotFoundError:
        print(f"File not found: {infile}")
    except Exception as e:
        print(f"An error occurred: {str(e)}")
    return taxids

def muscle_align_N_seq(seq_records: Iterator[SeqRecord]) -> Iterator[ SeqRecord ]:
    """Given a MSA of a protein class, and a fasta string of a chain, return a new MSA with the chain added to the class MSA."""
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_file:
        temp_filename = temp_file.name
        SeqIO.write([*seq_records], temp_filename, "fasta")
        muscle_cmd = [MUSCLE_BIN,'-in',temp_filename,'-quiet']

        try:
            process = subprocess.run(muscle_cmd, stdout=subprocess.PIPE, text=True)
            if  process.returncode == 0:
                muscle_out = process.stdout
                muscle_output_handle = StringIO(muscle_out)
                seq_records = SeqIO.parse(muscle_output_handle, "fasta")
                temp_file.close()
                os.remove(temp_filename)
                return seq_records
            else:
                print("{} failed with code {}".format(" ".join(muscle_cmd)), process.returncode)
                raise Exception("`{}` did not succeed.".format(" ".join(muscle_cmd)))
        except:
            temp_file.close()
            os.remove(temp_filename)
            raise Exception("Error running muscle.")


#! -------

#! -------
rib      = RibosomeAssets('3J7Z').profile()
organism_taxid = rib.src_organism_ids[0]
prots    = RibosomeAssets('3J7Z').profile().proteins

hmm_cachedir = ASSETS['__hmm_cache']

for candidate_class in list_ProteinClass:
    fasta_path   = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"{candidate_class}.fasta")
    records      = Fasta(fasta_path)
    ids          = records.all_taxids()
    phylo_nbhd   = phylogenetic_neighborhood(list(map(lambda x: str(x),ids)), str(organism_taxid), n_neighbors=10)
    seqs         = records.pick_taxids(phylo_nbhd)
    print("Picked seqnences,")
    print(seqs)
    seqs_aligned = muscle_align_N_seq( iter(seqs))
    seqs_aligned1, seqs_aligned2 = tee(seqs_aligned)    
    print("\n")
    print("Processing class {}".format(candidate_class))
    pprint([*seqs_aligned2])


    seq_tuples =  [TextSequence(name=bytes(seq.id, 'utf-8'), sequence=str(seq.seq)) for seq in seqs_aligned1]



    cached_name = "class_{}_taxid_{}.hmm".format(candidate_class, organism_taxid)

    alphabet      = pyhmmer.easel.Alphabet.amino() #* <---- AA?
    # alphabet      = pyhmmer.easel.Alphabet.rna() #* <---- NC?
    builder       = pyhmmer.plan7.Builder(alphabet)
    background    = pyhmmer.plan7.Background(alphabet) #? The null(background) model can be later augmented.
    anonymous_msa = pyhmmer.easel.TextMSA(bytes(cached_name, 'utf-8'),sequences=seq_tuples)

    hmm, _profile, _optmized_profile = builder.build_msa(anonymous_msa.digitize(alphabet), background)
    
    if not os.path.isfile(os.path.join(hmm_cachedir, cached_name)):
        with open(os.path.join(hmm_cachedir, cached_name), "wb") as hmm_file:
            hmm.write(hmm_file)
            print("Wrote `{}` to `{}`".format(cached_name, hmm_cachedir))
    else:
        print(os.path.join(hmm_cachedir, cached_name) + " exists. ")
    

# fasta_path = os.path.join(ASSETS["fasta_ribosomal_proteins"], f"bL9.fasta")
# record = [*SeqIO.parse(fasta_path, "fasta")]



#? extract sequences 
#? store (create an asset class in __init__), possibly compress
#? create method to pick phyl. nbhd, create an MSA, create an HMM
#? ->evaluate as before