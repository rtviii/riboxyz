"""We'd like to compare the conservation of a given region in the MSA against the background of a superkingdom's (or finer rank's) MSA"""
"""Protocol

I. Compare within superkingdom:
 - Get all available sequences (interpro): S
 - Reduce the same-species strains and sequences to a single representative consensus sequence
    Assumptions:
        - if the same species has fewer than 3 strains -- pick the longest
        - otherwise use dumb consensus (0.75 threshold)
 - align the species consesus seqs into superkingdom MSA:  K
 - construct an HMM of species' consensus sequences for the superkingdom 
 - asses the region of interest (e.g. constriction site) against the superkingdom HMM for each sequence:
    - { find outliers, insertions, deletions, reason from HMM stats... }



II. Compare between superkingdoms:
... WIP



|---------------- * ---------------- *---------------- * ---------------- Notes  ---------------- * ---------------- * ---------------- * ----------------|
* some of the species are still quite similar and therefore are overrepresented even after collapsing strains into species:

    e.g. Rotaria 
     2762511: {'members': [2762511, 2762511, 2762511, 2762511, 2762511, 2762511, 2762511, 2762511, 2762511, 2762511], 'name': 'Rotaria sp. Silwood1'},
     2762512: {'members': [2762512, 2762512, 2762512, 2762512, 2762512, 2762512, 2762512, 2762512, 2762512, 2762512, 2762512], 'name': 'Rotaria sp. Silwood2'},

  Perhaps some evolutionary metric can be used to "normalize" the tree over not just species<->strains relationships but along all ranks.
"""
# ?? A good first sanity check: find g.lambdia and the other 2 empirically found species via the hmm method.

from pprint import pprint
from typing import NewType
from ribctl.lib.libmsa import Fasta, Taxid, generate_consensus, ncbi, muscle_align_N_seq
from ribctl.lib import libhmm
from ribctl.lib import libmsa
from Bio.SeqRecord import SeqRecord

uL4_euk = Fasta("/home/rtviii/dev/riboxyz/__scripts/uL4_euk.fasta")

#! Consensus for substrains
def consensus_for_substrain(taxid):
    """This is poorly done. Both substrains and same-speices sequences get picked up """
    strains       = uL4_euk.pick_descendants_of_taxid(taxid)
    if len(strains) < 3:
        return (max(strains, key=lambda x: len(x.seq)))
    else:
        records       = muscle_align_N_seq(strains)
        consensus_seq = generate_consensus([*records]).dumb_consensus(threshold=0.75,require_multiple=True)
        return SeqRecord(seq=consensus_seq, description="Consensus sequence for taxids {}".format(" ".join(map(lambda r: r.id, strains)), id=taxid))


euk_species = uL4_euk.pick_descendants_of_taxid(libmsa.TAXID_EUKARYOTA)


def spec_strain_tree(_:list[SeqRecord])->dict:
    tree = {}
    for record in _:
        record_taxid  = int( record.id )
        species_taxid = Taxid.ancestor_at_rank(record_taxid, "species")

        if species_taxid not in tree:
            tree[species_taxid] = {"name": Taxid.get_name(species_taxid)[species_taxid], "members": [record]}
        else:
            tree[species_taxid]['members'].append(record)
    return tree

def normalize_from_tree(tree:dict)->list[SeqRecord]:
    """See notes above on full rank-normalization."""
    normalized  = []

    for species_taxid, v in tree.items():
        members = v['members']
        if len(members) > 2:
            aligned_  = muscle_align_N_seq(members)
            consensus = generate_consensus([*aligned_]).dumb_consensus(threshold=0.75,require_multiple=True, ambiguous='')
            normalized.append(SeqRecord(seq=consensus, description="Consensus sequence for taxids {}".format(" ".join(map(lambda r: r.id, v['members'])) ),id=str( species_taxid )))
        else:
            normalized.append(max(v['members'], key=lambda x: len(x.seq)))
    return normalized



norm = normalize_from_tree(spec_strain_tree(euk_species))

# Fasta(records = norm)
Fasta.write_fasta(norm, 'uL4_euk_normalized_tree.fasta')
pprint(norm)
pprint(len(norm))












exit()

# taxids = uL4_euk.all_taxids()
# euk, prok, arch = [], [], []
# eukaryotic_tree = {}

# # ! Partition into superkingdoms
# for t in taxidsrecord_taxidrecord_taxid:
#     supk = Taxid.superkingdom(t)
#     match (supk):
#         case "archaea":
#             arch.append(t)
#         case "eukaryota":
#             euk.append(t)
#         case "bacteria":
#             prok.append(t)

# euk_species = []

# for taxid in euk:
#     spec_id = Taxid.ancestor_at_rank(taxid, "species")
#     if spec_id not in eukarytic_tree:
#         eukaryotic_tree[spec_id] = {"name": Taxid.get_name(spec_id)[spec_id], "members": [taxid]}
#     else:
#         eukaryotic_tree[spec_id]['members'].append(taxid)
#         eukaryotic_tree[spec_id]['members'] = list(set(eukaryotic_tree[spec_id]['members']))


# for spec_id,v in eukaryotic_tree.items():
#     if len(v['members']) > 1:
#         print(spec_id, v['name'], v['members'])




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
