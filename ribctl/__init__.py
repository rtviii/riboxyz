import os
import pathlib
from typing import Literal


RIBETL_DATA = os.environ.get("RIBETL_DATA")


# This amounts to "_assets folder is expected to exist in the root of `ribctl`(next to top-level __init__.py)"
ASSETS_PATH = os.path.join(pathlib.Path(__file__).parent, "assets")

MUSCLE_BIN = os.path.join(ASSETS_PATH, "muscle3.8.1")
asset_type = Literal[
    "subunit_map_lsu",
    "subunit_map_ssu",
    "old_names_lsu",
    "old_names_ssu",

    "fasta_proteins_cytosolic",
    "fasta_proteins_mitochondrial",
    "fasta_ribosomal_rna",
    "fasta_factors",
    "fasta_trna",

    "__hmm_cache",
]

if os.environ.get("RIBETL_DATA") == "" or not os.path.exists(ASSETS_PATH):
    raise KeyError(
        "Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable."
    )

ASSETS: dict[asset_type, pathlib.Path] = {

    "subunit_map_lsu"         : pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_LSU.json")),
    "subunit_map_ssu"         : pathlib.Path(os.path.join(ASSETS_PATH, "subunit_map_SSU.json")),
    "old_names_lsu"           : pathlib.Path(os.path.join(ASSETS_PATH, "old_names_LSU.json")),
    "old_names_ssu"           : pathlib.Path(os.path.join(ASSETS_PATH, "old_names_SSU.json")),

    "fasta_proteins_cytosolic"    : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_ribosomal_proteins") ),
    "fasta_proteins_mitochondrial": pathlib.Path(os.path.join(ASSETS_PATH, "fasta_ribosomal_proteins") ),
    "fasta_ribosomal_rna"         : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_ribosomal_rna") ),
    "fasta_factors"               : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_factors")),
    "fasta_trna"                  : pathlib.Path(os.path.join(ASSETS_PATH, "fasta_trna")),

    "__hmm_cache"             : pathlib.Path(os.path.join(ASSETS_PATH, "__hmm_cache")),
}

TAXID_BACTERIA          = 2
TAXID_EUKARYA           = 2759
TAXID_ARCHEA            = 2157

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

model_species = {
    "Lactococcus"               : 1357,
    "Mycobacterium smegmatis"   : 1772,
    "Mycobacterium tuberculosis": 1773,
    "Bacillus subtilis"         : 1423,
    "Leishmania donovani"       : 5661,
    "Trypanosoma cruzi"         : 5693,
    "Trichomonas vaginalis"     : 5722,
    "Giardia duodenalis"        : 5741,
    "Spraguea lophii"           : 51541,
}

model_subgenuses = {
    "Lactococcus"               : 1357,
    "Mycobacterium smegmatis"   : 1772,
    "Mycobacterium tuberculosis": 1773,
    "Bacillus subtilis"         : 1423,
    "Leishmania donovani"       : 38568,
    "Trypanosoma cruzi"         : 47570,
    "Trichomonas vaginalis"     : 181550,
    "Giardia duodenalis"        : 5741,
    "Spraguea lophii"           : 51540,
}


__version__ = "0.1.0"
