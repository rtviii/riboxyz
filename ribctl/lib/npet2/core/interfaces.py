# ribctl/lib/npet2/core/interfaces.py
"""
Provider protocols for npet2.

These define the boundary between the pipeline and any data source.
Implement these to plug in riboxyz, a local file system, or an API.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Protocol

import numpy as np

from .ribosome_types import RibosomeProfile, PTCInfo, ConstrictionInfo
from .types import ArtifactRef, ArtifactType


class StructureProvider(Protocol):
    """Provides atom coordinates and the ribosome profile for a structure."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        """
        Must return:
          - atom_xyz: (N, 3) float32
          - atom_element: (N,) str array (optional)
          - mmcif_path: str
          - profile: RibosomeProfile
        """
        ...


class LandmarkProvider(Protocol):
    """Provides PTC and constriction site coordinates."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        """
        Must return:
          - ptc_xyz: (3,) float32
          - constriction_xyz: (3,) float32
        """
        ...


class ArtifactStore(Protocol):
    @property
    def run_dir(self) -> Path: ...

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType,
                  data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_json(self, *, name: str, stage: str, obj: Any,
                 meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray,
                  meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def add_ref(self, ref: ArtifactRef) -> None: ...

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None: ...
