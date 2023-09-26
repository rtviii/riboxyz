from Bio import SeqIO
from ete3 import NCBITaxa

# Initialize the NCBI Taxonomy database
ncbi = NCBITaxa()


# Translate the taxonomic name to a taxonomic ID
# Define the taxonomic name you want to translate

# Replace 'your_fasta_file.fasta' with the path to your FASTA file
fasta_file = "./16S_mitochondrial_rRNA_AND_so_rna_type_nameMt_rRNA_AND_entry_typeSequence.fasta"
dest       = "./m_16SrRNA.fasta"

to_keep    = []


# Use SeqIO to read the input FASTA file and store records in the list
# Save the records in a new FASTA file

# Use SeqIO to open and iterate through the FASTA file
i = 0
specs = {}

with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        # if "16S" in record.description and ( len(record.seq) < 1200 and len(record.seq) > 600 ):
            try:
                taxname_huh = " ".join(record.description.split(" ")[1:3])
                taxid       = str(ncbi.get_name_translator([taxname_huh])[taxname_huh][0])
                record.id = str(taxid)

                if taxid not in specs:
                    specs[taxid] = 1
                else:
                    specs[taxid] += 1

                record.description = "rnacentral_id:{}|{}".format(record.id, record.description)

                if specs[taxid] < 3:
                    to_keep.append(record)
            except:
                ...

with open(dest, "w") as output_handle:
    SeqIO.write(to_keep, output_handle, "fasta")