# ribctl/lib/npet2/adapters/standalone_providers.py
"""
Standalone providers for running npet2 without the riboxyz repo.

Input sources:
  - mmCIF file on disk
  - Profile JSON (from riboxyz API or a local file)
  - Landmark coordinates (from riboxyz API or a local file)
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from Bio.PDB.MMCIFParser import FastMMCIFParser

from ribctl.lib.npet2.core.ribosome_types import (
    ConstrictionInfo,
    PTCInfo,
    RibosomeProfile,
)
from ribctl.lib.npet2.core.settings import RIBOXYZ_API_BASE


class FileStructureProvider:
    """
    Load atoms from a local mmCIF file and profile from a JSON file or API.

    Usage:
        provider = FileStructureProvider(
            mmcif_path="/data/7K00.cif",
            profile_path="/data/7K00_profile.json",  # or None to fetch from API
        )
    """

    def __init__(
        self,
        mmcif_path: str | Path,
        profile_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.mmcif_path = Path(mmcif_path)
        self.profile_path = Path(profile_path) if profile_path else None
        self.api_base = api_base or RIBOXYZ_API_BASE

        if not self.mmcif_path.exists():
            raise FileNotFoundError(f"mmCIF file not found: {self.mmcif_path}")

    def fingerprint(self, rcsb_id: str) -> str:
        return f"mmcif:{self.mmcif_path}"

    def _load_profile(self, rcsb_id: str) -> RibosomeProfile:
        if self.profile_path and self.profile_path.exists():
            data = json.loads(self.profile_path.read_text())
            return RibosomeProfile.model_validate(data)

        # Fetch from API
        import requests

        url = f"{self.api_base}/structures/{rcsb_id.upper()}/profile"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return RibosomeProfile.model_validate(resp.json())

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        parser = FastMMCIFParser(QUIET=True)
        structure = parser.get_structure(rcsb_id, str(self.mmcif_path))
        atoms = list(structure[0].get_atoms())

        if not atoms:
            raise ValueError(f"No atoms found in {self.mmcif_path}")

        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray(
            [getattr(a, "element", "") or a.get_id()[0] for a in atoms]
        )

        profile = self._load_profile(rcsb_id)

        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(self.mmcif_path),
            "profile": profile,
        }


class FileLandmarkProvider:
    """
    Load PTC/constriction from a local JSON file or the riboxyz API.

    Local file format:
        {
            "ptc": {"location": [x, y, z]},
            "constriction": {"location": [x, y, z]}
        }
    """

    def __init__(
        self,
        landmarks_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.landmarks_path = Path(landmarks_path) if landmarks_path else None
        self.api_base = api_base or RIBOXYZ_API_BASE

    def fingerprint(self, rcsb_id: str) -> str:
        if self.landmarks_path:
            return f"landmarks_file:{self.landmarks_path}"
        return f"landmarks_api:{self.api_base}"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        if self.landmarks_path and self.landmarks_path.exists():
            data = json.loads(self.landmarks_path.read_text())
            ptc_info = PTCInfo.model_validate(data["ptc"])
            constr_info = ConstrictionInfo.model_validate(data["constriction"])
            return {
                "ptc_xyz": np.array(ptc_info.location, dtype=np.float32),
                "constriction_xyz": np.array(constr_info.location, dtype=np.float32),
            }

        # Fetch from API
        import requests

        rcsb_id = rcsb_id.upper()

        ptc_resp = requests.get(
            f"{self.api_base}/loci/ptc", params={"rcsb_id": rcsb_id}, timeout=30
        )
        ptc_resp.raise_for_status()
        ptc_info = PTCInfo.model_validate(ptc_resp.json())

        constr_resp = requests.get(
            f"{self.api_base}/loci/constriction_site",
            params={"rcsb_id": rcsb_id},
            timeout=30,
        )
        constr_resp.raise_for_status()
        constr_info = ConstrictionInfo.model_validate(constr_resp.json())

        return {
            "ptc_xyz": np.array(ptc_info.location, dtype=np.float32),
            "constriction_xyz": np.array(constr_info.location, dtype=np.float32),
        }