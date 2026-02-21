# ribctl/lib/exceptions.py

class SkipAsset(RuntimeError):
    """
    Raised when an asset is not applicable for a given structure.

    Examples:
      - SSU-only entry: no uL4/uL22 -> constriction not defined
      - SSU-only entry: no LSU rRNA -> PTC not defined
      - Assembly intermediates: LSU not present/annotated -> skip
      - Mapping fails due to missing reference residues -> skip (not a hard failure)
    """
    pass
