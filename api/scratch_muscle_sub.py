from prody import MSA, calcShannonEntropy
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import RNA, PolymericFactor, Protein, ProteinClass
from api.ribctl.msa.msalib import msa_class_proteovision_path, msa_profiles_dict, msa_profiles_dict_prd, prot_class_msa_extend, prot_class_msa_extend_prd





rcsb_id    :str         = "5AFI"
poly_class:ProteinClass = "bS6"
R = RibosomeAssets(rcsb_id)
chain             = R.get_chain_by_polymer_class(poly_class)
if chain is None:raise LookupError() 

new_seq  = chain.entity_poly_seq_one_letter_code_can.replace("\n","")
msa_profiles: dict[ProteinClass, MSA] = msa_profiles_dict_prd()

for class_name, msa in list( msa_profiles.items() ):
    extended_class = prot_class_msa_extend_prd(rcsb_id,class_name, new_seq)
    H_extended     = sum(calcShannonEntropy(extended_class))
    H_original     = sum(calcShannonEntropy(msa))
    H_delta        = H_extended - H_original
    print(class_name, "\t:", H_delta)


# TODO: Add an phylogeny adjustment if organism is provided.
async def seq_asses_protclass_fit(sequence:str,prot_class_profile:MSA )->float:

    extended_class = prot_class_msa_extend_prd(rcsb_id,class_name, new_seq)
    H_extended     = sum(calcShannonEntropy(extended_class))
    H_original     = sum(calcShannonEntropy(msa))
    H_delta        = H_extended - H_original
    print(class_name, "\t:", H_delta)
    return 0
    