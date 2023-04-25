import os
from pprint import pprint
import subprocess
from Bio import SeqRecord
from prody import Polymer
from Bio.Seq import Seq
from Bio import AlignIO
from io import StringIO
import prody
from api.ribctl.etl.ribosome_assets import RibosomeAssets
from api.ribctl.lib.types.types_ribosome import RNA, PolymericFactor, Protein, ProteinClass
from api.ribctl.msa.msalib import msa_class_proteovision_path, prot_class_msa_extend, seq_to_fasta



rcsb_id    :str         = "5AFI"
poly_class:ProteinClass = "bL25"


msa = prot_class_msa_extend(rcsb_id,poly_class)
