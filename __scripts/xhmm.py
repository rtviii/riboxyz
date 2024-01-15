"""Another one is create an MSA for a given rank (e.g. superkingdom) and compare the conservation of a given region in the MSA against this background"""
"""Protocol
Within superkingdom:
 - Get all available sequences (interpro): S
 - Reduce the same-species strains to a single representative consensus seq 
 - align the species consesus seqs into superkingdom MSA:  K
 - construct an HMM of species' consensus sequences for the superkingdom 
 - asses the region of interest (e.g. constriction site) against the superkingdom HMM for each sequence:
    - { find outliers, insertions, deletions, reason from HMM stats... }

Compare between superkingdoms:
...

"""
# ?? A good first sanity check: find g.lambdia and the other 2 empirically found species.

from pprint import pprint
from ribctl.lib.libmsa import Fasta, Taxid, ncbi

uL4_euk = Fasta("/home/rtviii/dev/riboxyz/__scripts/uL4_euk.fasta")

# f = Fasta("/home/rtviii/dev/riboxyz/ribctl/uL22-protein-matching-IPR001063.fasta")

taxids = uL4_euk.all_taxids()
euk, prok, arch = [], [], []
eukaryotic_tree = {}

# ! Partition into superkingdoms
for t in taxids:
    supk = Taxid.superkingdom(t)
    match (supk):
        case "archaea":
            arch.append(t)
        case "eukaryota":
            euk.append(t)
        case "bacteria":
            prok.append(t)

euk_species = []

for taxid in euk:
    spec_id = Taxid.ancestor_at_rank(taxid, "species")
    if spec_id not in eukaryotic_tree:
        eukaryotic_tree[spec_id] = {"name": Taxid.get_name(spec_id)[spec_id], "members": [taxid]}
    else:
        eukaryotic_tree[spec_id]['members'].append(taxid)
        eukaryotic_tree[spec_id]['members'] = list(set(eukaryotic_tree[spec_id]['members']))


for spec_id,v in eukaryotic_tree.items():
    if len(v['members']) > 1:
        print(spec_id, v['name'], v['members'])



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
