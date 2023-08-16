import os

RIBETL_DATA = os.environ.get("RIBETL_DATA") 

if os.environ.get("RIBETL_DATA") == "":
    raise KeyError("Repostiry of static PDB files should be defined as $RIBETL_DATA environment variable.")


# __all__ = [...]

__version__ = '0.1.0'
