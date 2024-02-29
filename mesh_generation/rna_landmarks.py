# The idea is to see if there are enough conserved residues in the LSU rRNA in the vicinity of the uL4/uL22 tips:
# grab N closest pairs of atoms from the uL4/uL22 tips and then look at the alignment of rRNA for this kingdom a
# see if enough conserved columns are in the vicinity (increasing radius) to form an alpha shape around the constriction site.



from pprint import pprint
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libhmm import fasta_phylogenetic_correction
from ribctl.lib.libmsa import Fasta, muscle_align_N_seq
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass


RCSB_ID = "4W29"
ra = RibosomeAssets(RCSB_ID)

c         = ra.get_chain_by_polymer_class("23SrRNA")
src_taxid = ra.get_taxids()[0][0]
seqs      = list(fasta_phylogenetic_correction(PolymerClass("23SrRNA"), src_taxid))

seqs_a = muscle_align_N_seq(seqs,vvv=True)
Fasta.write_fasta(list(seqs_a), "23SrRNA_neighbors_of_{}.fasta".format(RCSB_ID))