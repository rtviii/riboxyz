import os
import pathlib
from typing import Literal

def get_env_var(name: str) -> str:
    value = os.environ.get(name)
    if not value:
        raise KeyError(f"CRITICAL: Environment variable '{name}' is not set in .env.")
    return value

# Mandatory Directory Checks
RIBETL_DATA = pathlib.Path(get_env_var("RIBETL_DATA"))
if not RIBETL_DATA.exists():
    raise NotADirectoryError(f"RIBETL_DATA path does not exist: {RIBETL_DATA}")
#! ------------- logs ----------------
CLASSIFICATION_REPORTS = os.path.join(pathlib.Path(__file__).parent, "logs","hmm_classification_reports")

#! ------------- assets ----------------
# This amounts to "assets folder is expected to exist in the root of `ribctl`(next to top-level __init__.py)"
ASSETS_PATH        = os.path.join(pathlib.Path(__file__).parent, "assets_project")
RIBXZ_TEMP_FILES   = os.path.join(ASSETS_PATH,"temp")
MUSCLE_BIN         = os.path.join(ASSETS_PATH, "muscle3.8.1")
NCBI_TAXDUMP_GZ    = os.path.join(ASSETS_PATH, "taxdump.tar.gz")
NCBI_TAXA_SQLITE = os.environ.get("NCBI_TAXA_SQLITE")


# Mandatory File Checks (No more fallbacks to ASSETS_PATH)
NCBI_TAXA_SQLITE = get_env_var("NCBI_TAXA_SQLITE")
if not os.path.exists(NCBI_TAXA_SQLITE):
    raise FileNotFoundError(f"NCBI_TAXA_SQLITE file not found at: {NCBI_TAXA_SQLITE}")

# Assets (If these are required for production, treat them the same)
ASSETS_PATH = get_env_var("ASSETS_PATH")
CHAINSPLITTER_PATH = os.path.join(pathlib.Path(__file__).parent.parent, "chimerax", "chainsplitter.py")

#! -------------- locations
EXIT_TUNNEL_WORK  = os.path.join(ASSETS_PATH, "exit_tunnel_work")
POISSON_RECON_BIN = os.path.join(EXIT_TUNNEL_WORK, "PoissonRecon")

asset_type  = Literal[
    "subunit_map_lsu",
    "subunit_map_ssu",
    "old_names_lsu",
    "old_names_ssu",

    "fasta_proteins_cytosolic",
    "fasta_proteins_mitochondrial",
    "fasta_rna",

    "fasta_factors_elongation",
    "fasta_factors_elongation_archaea",
    "fasta_factors_elongation_bacteria",
    "fasta_factors_elongation_eukaryota",

    "fasta_factors_initiation",
    "fasta_factors_initiation_archaea",
    "fasta_factors_initiation_bacteria",
    "fasta_factors_initiation_eukaryota",

    "cache_hmm",
]

# if os.environ.get("RIBETL_DATA") == "" or not os.path.exists(ASSETS_PATH):
#     raise KeyError(
#         "Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable.")
if os.environ.get(NCBI_TAXA_SQLITE) == "" or not os.path.exists(NCBI_TAXA_SQLITE):
    import warnings
    warnings.warn("""NCBI taxonomy sqlite file should be available at NCBI_TAXA_SQLITE environment variable. 
        The dump will be downloaded and unpacked by ete3 automatically (from https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz).
                   """)


# if os.envion.get(NCBI_TAXDUMP_GZ) == "" or not os.path.exists(NCBI_TAXDUMP_GZ):
#     raise FileNotFoundError(
#         """NCBI taxonomy dump file should be available at $NCBI_TAXDUMP_GZ environment variable. 
#         Download it here: https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz """
#     )

ASSETS: dict[asset_type, pathlib.Path] = {
    "subunit_map_lsu"         : pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_LSU.json")),
    "subunit_map_ssu"         : pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_SSU.json")),
    "old_names_lsu"           : pathlib.Path(os.path.join(ASSETS_PATH, "old_names_LSU.json")),
    "old_names_ssu"           : pathlib.Path(os.path.join(ASSETS_PATH, "old_names_SSU.json")),

    "fasta_proteins_cytosolic"    : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_proteins_cytosolic") ),
    "fasta_proteins_mitochondrial": pathlib.Path(os.path.join(ASSETS_PATH, "fasta_proteins_mitochondrial") ),

    "fasta_rna": pathlib.Path(os.path.join(ASSETS_PATH, "fasta_rna") ),

    "fasta_factors_initiation"    : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_initiation")),
    "fasta_factors_initiation_archaea"  : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_initiation","archaea")),
    "fasta_factors_initiation_bacteria" : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_initiation","bacteria")),
    "fasta_factors_initiation_eukaryota": pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_initiation","eukaryota")),

    "fasta_factors_elongation"          : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_elongation")),
    "fasta_factors_elongation_bacteria" : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_elongation","bacteria")),
    "fasta_factors_elongation_eukaryota": pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_elongation","eukaryota")),
    "fasta_factors_elongation_archaea"  : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors_elongation","archaea")),

    "cache_hmm"             : pathlib.Path(os.path.join(ASSETS_PATH, "cache_hmm")),
}

TAXID_BACTERIA  = 2
TAXID_EUKARYOTA = 2759
TAXID_ARCHAEA   = 2157

AMINO_ACIDS_3_TO_1_CODE = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "ASX": "B",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLX": "Z",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
AMINO_ACIDS_1_TO_3_CODE = {v: k for k, v in AMINO_ACIDS_3_TO_1_CODE.items()}

__version__ = "0.1.0"
