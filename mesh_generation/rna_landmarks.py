# The idea is to see if there are enough conserved residues in the LSU rRNA in the vicinity of the uL4/uL22 tips:
# grab N closest pairs of atoms from the uL4/uL22 tips and then look at the alignment of rRNA for this kingdom a
# see if enough conserved columns are in the vicinity (increasing radius) to form an alpha shape around the constriction site.



from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libmsa import muscle_align_N_seq


RCSB_ID = "4W29"
ra = RibosomeAssets(RCSB_ID)

c = ra.get_chain_by_polymer_class("23SrRNA")
print(c)

# muscle_align_N_seq()