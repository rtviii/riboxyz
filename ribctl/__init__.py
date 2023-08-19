import os
import pathlib
from typing import Literal


RIBETL_DATA = os.environ.get("RIBETL_DATA")

# This amounts to "_assets folder is expected to exist in the root of `ribctl`(next to top-level __init__.py)"
_ASSETS_PATH  = os.path.join(pathlib.Path(__file__).parent, '_assets')
asset_type =  Literal[
                    "landmark_sites",
                    "hmm_ribosomal_proteins",
                    "hmm_factors",
                    "hmm_ribosomal_rna",
                    "hmm_trna"
                    ]

if os.environ.get("RIBETL_DATA") == "" or not os.path.exists(_ASSETS_PATH):
    raise KeyError("Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable.")

ASSETS:dict[asset_type, pathlib.Path] = {
'landmark_sites'        : pathlib.Path(os.path.join(_ASSETS_PATH, 'landmark_sites')),
'hmm_ribosomal_proteins': pathlib.Path(os.path.join(_ASSETS_PATH, 'hmm_ribosomal_proteins')),
'hmm_factors'           : pathlib.Path(os.path.join(_ASSETS_PATH, 'hmm_factors')),
'hmm_ribosomal_rna'     : pathlib.Path(os.path.join(_ASSETS_PATH, 'hmm_ribosomal_rna')),
'hmm_trna'              : pathlib.Path(os.path.join(_ASSETS_PATH, 'hmm_trna'))
}



__version__ = '0.1.0'
