# The idea is to see if there are enough conserved residues in the LSU rRNA in the vicinity of the uL4/uL22 tips:
# grab N closest pairs of atoms from the uL4/uL22 tips and then look at the alignment of rRNA for this kingdom a
# see if enough conserved columns are in the vicinity (increasing radius) to form an alpha shape around the constriction site.

import json
from pprint import pprint
import numpy as np
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libhmm import fasta_phylogenetic_correction
from ribctl.lib.libmsa import Fasta, muscle_align_N_seq
from ribctl.lib.ribosome_types.types_ribosome import PolymerClass
from Bio.PDB.Chain import Chain
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import Selection
from scipy.spatial.distance import cdist

import urllib.request as urlreq
from dash import Dash, html
import dash_bio as dashbio
from paths import tunnel_atom_encoding_path


RCSB_ID = "4W29"
ra = RibosomeAssets(RCSB_ID)

# # !
# c         = ra.get_chain_by_polymer_class("23SrRNA")
# seq_c = c.to_SeqRecord()
# src_taxid = ra.get_taxids()[0][0]
# seqs      = list(fasta_phylogenetic_correction(PolymerClass("23SrRNA"), src_taxid))

# seqs_a = muscle_align_N_seq([ seq_c, *seqs ],vvv=True)
# Fasta.write_fasta(list(seqs_a), "23SrRNA_neighbors_of_{}.fasta".format(RCSB_ID))

FASTA_PATH  = '/home/rtviii/dev/riboxyz/23SrRNA_neighbors_of_4W29.fasta'
msa          = Fasta(FASTA_PATH)


from ribctl.lib.libmsa import util__backwards_match, util__forwards_match

with open(tunnel_atom_encoding_path(RCSB_ID), 'r') as infile:
    encoding_data = json.load(infile)

rna_residues =[]
for atom in encoding_data:
    if len(atom['chain_nomenclature'] ) > 0:
        if 'RNA' in atom['chain_nomenclature'][0]:
            rna_residues.append(atom)

vicinity_residues_seqids = sorted(list(set([a['residue_seqid'] for a  in rna_residues])))

rna_record_a = None

for record in msa.records:
    if RCSB_ID in record.description:
       rna_record_a = record

aligned_column_ix = []
for res_seqid in vicinity_residues_seqids:
    ix = util__forwards_match(rna_record_a._seq, res_seqid)
    aligned_column_ix.append(ix)

def nucleotides_at_columns(msa: Fasta, columns_indices: list[int]) -> list[list[str]]:
    _ = []
    for index in columns_indices:
        col_ =[]
        for record in msa.records:
                col_.append(record[index])
        _.append([index,col_])
        
    return _


pprint(nucleotides_at_columns(msa, aligned_column_ix))
pprint(len(nucleotides_at_columns(msa, aligned_column_ix)))



# pprint(aligned_column_ix)
# pprint(vicinity_residues_seqids)

#TODO: "Get column" from the alignment (msa)
#TODO: visualize conserved
#TODO: calculate conservation for every column
#TODO: assign distance from centerline to every column
#TODO: HMM shit


# looking at the structure

if False: # VIA CONSTRICTION SITE
    ul4     = ra.get_chain_by_polymer_class('uL4')
    ul22    = ra.get_chain_by_polymer_class('uL22')
    LSUrRNA = ra.get_LSU_rRNA()


    if None in [ul22, ul4, LSUrRNA]:
        raise Exception("One of the chains is None: ",  [ul22, ul4, LSUrRNA])

    structure = ra.biopython_structure()

    for model in structure:
        for chain in model:
            # Check if the chain's auth_asym_id matches the provided value
            if chain.id == ul4.auth_asym_id:
                ul4_chain:Chain = chain
            if chain.id == ul22.auth_asym_id:
                ul22_chain:Chain = chain
            if chain.id == LSUrRNA.auth_asym_id:
                LSUrRNA_chain:Chain = chain


    coordinates_uL4  = np.array([np.array(atom.get_coord()) for atom in ul4_chain.get_atoms()])
    coordinates_uL22 = np.array([np.array(atom.get_coord()) for atom in ul22_chain.get_atoms()])

    distances = cdist(coordinates_uL4, coordinates_uL22, 'euclidean')
    min_index = np.unravel_index(np.argmin(distances), distances.shape)
    closest_point_cloud1 = coordinates_uL4[min_index[0]]
    closest_point_cloud2 = coordinates_uL22[min_index[1]]
    constriction_centroid_uL22_uL4 = np.mean([closest_point_cloud1,closest_point_cloud2], axis=0)

    #! further:
    # grab the 23S rna atoms in the vicinity of the centroid
    # rrna_residues = [*LSUrRNA_chain.get_residues()]
    rrna_residues = [*LSUrRNA_chain.get_atoms()]
    print(rrna_residues)
    # NOw locate the ones within 5/10/15 A of the centroid

    centroid_nbhd_atoms = []
    ns            = NeighborSearch(list(structure.get_atoms()))
    nearby_atoms  = ns.search(constriction_centroid_uL22_uL4,10, "A")
    centroid_nbhd_atoms.extend(nearby_atoms)

