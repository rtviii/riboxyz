
# muscle_combine_profiles(chief_msa_path, control_path, '{}with_{}.fasta'.format(chief_msa_path,control))
# h, hwctl = compare_entropies(chief_msa_path, '{}with_{}.fasta'.format(chief_msa_path,control))

# ※ ---------------------
# Goals are to help Shiqi find the the common landmark for the exit port:
    # For every domain:
# 1. Pick the target in one structure X (let's say it's (nucleotide U61 in chain ??) or (amino-acid ?? in uL23)). Call this the origin site: chain + nucleotide/amino-acid.
# 2. Given the species of X, get an alignment of all the other species in this domain. 
# 3. Align the original chain into this MSA (profile-profile in muscle). Call the resulting [columns] aligned to the origin site in the MSA: landmark columns/indices
        # For all other structures of interest in this domain:
        # 4. Get the chain of the same class and align it into the MSA.
        # 5. Track back the reference columns from the MSA-stretched chain to the original (count and exclude the gaps basically)