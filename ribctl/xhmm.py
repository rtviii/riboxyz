"""Another one is create an MSA for a given rank (e.g. superkingdom) and compare the conservation of a given region in the MSA against this background"""

from ribctl.lib.msalib import Fasta, Taxid

f = Fasta("/home/rtviii/dev/riboxyz/ribctl/uL22-protein-matching-IPR001063.fasta")


taxids = f.all_taxids(lambda x: int( (x.split("|")[-1]).split(":")[-1]) )
euk, prok, arch = [], [], []

for t in taxids:
    supk = Taxid.superkingdom(t)
    print(supk)
    match (supk):
        case 'archaea':
            arch.append(t)
        case 'eukaryota':
            euk.append(t)
        case 'bacteria':
            prok.append(t)



print(len(euk), len(prok), len(arch))

# Fasta.fasta_display_species(taxids)





# """One simple approach is just to grab the conserved region and based on the shannon entropy determine the deviation of a given
# member of the MSA"""
# from Bio import AlignIO
# import numpy as np
# from scipy import stats
# import matplotlib.pyplot as plt

# # Step 1: Read FASTA file and perform alignment (if necessary)
# alignment = AlignIO.read("your_file.fasta", "fasta")

# # Step 2: Extract region of interest (200-250)
# region_of_interest = alignment[:, 199:250]

# # Step 3: Calculate conservation scores (using Shannon entropy as an example)
# def shannon_entropy(column):
#     _, counts = np.unique(column, return_counts=True)
#     probabilities = counts / len(column)
#     return -np.sum(probabilities * np.log2(probabilities))

# conservation_scores = [shannon_entropy(column) for column in region_of_interest.T]

# # Step 4: Perform statistical tests (e.g., Chi-square test)
# # You'll need to define expected frequencies based on background or reference data.

# # Step 5: Visualization
# plt.plot(range(200, 251), conservation_scores)
# plt.xlabel('Position')
# plt.ylabel('Conservation Score')
# plt.title('Conservation Scores for Region of Interest')
# plt.show()
