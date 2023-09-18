import os
import pathlib
from typing import Literal


RIBETL_DATA = os.environ.get("RIBETL_DATA")


# This amounts to "_assets folder is expected to exist in the root of `ribctl`(next to top-level __init__.py)"
ASSETS_PATH  = os.path.join(pathlib.Path(__file__).parent, 'assets')

#TODO: make msas bona-fide assets
RP_MSAS_PATH        = os.path.join(ASSETS_PATH, 'rp_class_msas')
RP_MSAS_PRUNED_PATH = os.path.join(ASSETS_PATH, 'rp_class_msas_pruned')

asset_type =  Literal[
                    "subunit_map_lsu",
                    "subunit_map_ssu",
                    "old_names_lsu",
                    "old_names_ssu",

                    "fasta_ribosomal_proteins",
                    "fasta_ribosomal_rna",
                    "fasta_factors",
                    "fasta_trna",
                    
                    "hmm_cache",
                    "hmm_ribosomal_proteins",
                    "hmm_factors",
                    "hmm_ribosomal_rna",
                    "hmm_trna",
                    ]

if os.environ.get("RIBETL_DATA") == "" or not os.path.exists(ASSETS_PATH):
    raise KeyError("Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable.")

ASSETS:dict[asset_type, pathlib.Path] = {
 
"subunit_map_lsu":  pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_LSU.json")),
"subunit_map_ssu":  pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_SSU.json")),
"old_names_lsu":  pathlib.Path(os.path.join(ASSETS_PATH, "old_names_LSU.json")),
"old_names_ssu":  pathlib.Path(os.path.join(ASSETS_PATH, "old_names_SSU.json")),


'fasta_ribosomal_proteins': pathlib.Path(os.path.join(ASSETS_PATH, "fasta_ribosomal_proteins")),
'fasta_ribosomal_rna'     : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_ribosomal_rna")),
'fasta_factors'           : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors")),
'fasta_trna'              : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_trna")),

'hmm_ribosomal_proteins'  : pathlib.Path(os.path.join(ASSETS_PATH, 'hmm_ribosomal_proteins')),
'hmm_factors'             : pathlib.Path(os.path.join(ASSETS_PATH, 'hmm_factors')),
'hmm_ribosomal_rna'       : pathlib.Path(os.path.join(ASSETS_PATH, 'hmm_ribosomal_rna')),
'hmm_trna'                : pathlib.Path(os.path.join(ASSETS_PATH, 'hmm_trna'))
}

TAXID_BACTERIA          = 2
TAXID_EUKARYA           = 2759
TAXID_ARCHEA            = 2157
AMINO_ACIDS_3_TO_1_CODE = {
"ALA":"A",
"ARG":"R",
"ASN":"N",
"ASP":"D",
"ASX":"B",
"CYS":"C",
"GLU":"E",
"GLN":"Q",
"GLX":"Z",
"GLY":"G",
"HIS":"H",
"ILE":"I",
"LEU":"L",
"LYS":"K",
"MET":"M",
"PHE":"F",
"PRO":"P",
"SER":"S",
"THR":"T",
"TRP":"W",
"TYR":"Y",
"VAL":"V"};
AMINO_ACIDS_1_TO_3_CODE = {v:k for k,v in AMINO_ACIDS_3_TO_1_CODE.items()}



__version__ = '0.1.0'
