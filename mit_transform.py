from Bio import SeqIO
from ete3 import NCBITaxa

# Initialize the NCBI Taxonomy database
ncbi = NCBITaxa()


# Translate the taxonomic name to a taxonomic ID
# Define the taxonomic name you want to translate

# Replace 'your_fasta_file.fasta' with the path to your FASTA file
fasta_file = "./12S_mitochondrial_rRNA_AND_entry_typeSequence.fasta"

# Use SeqIO to open and iterate through the FASTA file
i = 0
with open(fasta_file, "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
        if "12S" in record.description:
            try:
                i=i+1
                taxname_huh = " ".join(record.description.split(" ")[1:3])
                taxid = ncbi.get_name_translator([taxname_huh])[taxname_huh]
                # print(taxid, len(record.seq))
                # print(record.description)
            except:
                ...
        else:
            print("\033[31m {}  \033[0m".format(record.description))


        # print("Sequence:", record.seq)
        # print("Length:", len(record))

        # Perform any additional processing or analysis here
print("Got 12S seqs : ", i)