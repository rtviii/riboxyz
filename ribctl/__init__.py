import os
import pathlib
from typing import Literal


RIBETL_DATA = os.environ.get("RIBETL_DATA")
ASSET_PATH  = os.path.join(pathlib.Path(__file__), '_assets')

_assets =  Literal["landmark_sites", "hmm_ribosomal_proteins", "hmm_factors", "hmm_ribosomal_rna", "hmm_trna"]





if os.environ.get("RIBETL_DATA") == "":
    raise KeyError("Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable.")


# __all__ = [...]

__version__ = '0.1.0'
