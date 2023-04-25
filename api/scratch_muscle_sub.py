from collections import defaultdict
from scipy.stats import entropy
import collections
from io import StringIO
from Bio import pairwise2
import os
from pprint import pprint
from prody import MSA, MSAFile, buildMSA
from Bio import AlignIO
from prody import calcShannonEntropy
import prody
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import RNA, PolymericFactor, Protein, ProteinClass
from api.ribctl.msa.msalib import msa_class_proteovision_path, msa_profiles_dict, msa_profiles_dict_prd, prot_class_msa_extend



# similarity_scores = defaultdict(float)




rcsb_id    :str         = "5AFI"
poly_class:ProteinClass = "bL25"
R = RibosomeAssets(rcsb_id)
chain             = R.get_chain_by_polymer_class(poly_class)
if chain is None:raise LookupError() 
new_seq = chain.entity_poly_seq_one_letter_code_can.replace("\n","")


# msa_dict          = msa_profiles_dict()
msa_dict          = msa_profiles_dict()
similarity_scores = defaultdict(float)


def column_entropy(seq:str):
    ntds = collections.Counter([ntd for ntd in seq])
    dist = [x/sum(ntds.values()) for x in ntds.values()]
    H = entropy(dist, base=2)
    return H

for class_name, msa in msa_dict.items():



   

        
    for i in len(list(msa)):
        column_entropy(msa[:,i])








    # extended = prot_class_msa_extend(rcsb_id,poly_class)
# 


    # print(class_name, H)
    # H_w_ctl      = sum(calcShannonEntropy(msa_file_ctl, omitgaps=False))

    # alignment = pairwise2.align.globalxx(new_seq, seq_record.seq)
    # similarity_score = alignment[0].score / len(alignment[0].seqA)
    # similarity_scores[class_name] += similarity_score

# pprint(similarity_scores)




# print(len(msa_profiles))
# print(msa_profiles)

# for class_name, msa in msa_dict.items():
#     for seq_record in msa:
#         alignment = pairwise2.align.globalxx(new_seq, seq_record.seq)
#         similarity_score = alignment[0].score / len(alignment[0].seqA)
#         similarity_scores[class_name] += similarity_score