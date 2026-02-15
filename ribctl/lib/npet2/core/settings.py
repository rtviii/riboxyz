# ribctl/lib/npet2/core/settings.py
"""
Default paths and external tool locations for npet2.

All of these can be overridden via:
  - Environment variables (NPET2_ROOT, NPET2_POISSON_RECON_BIN, etc.)
  - CLI flags
  - RunConfig / programmatic construction
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(var: str, default: str) -> Path:
    return Path(os.environ.get(var, default))


NPET2_ROOT = _env_path("NPET2_ROOT", str(Path.home() / "npet2_data"))
NPET2_RUNS_ROOT = _env_path("NPET2_RUNS_ROOT", str(NPET2_ROOT / "runs"))
NPET2_CACHE_ROOT = _env_path("NPET2_CACHE_ROOT", str(NPET2_ROOT / "cache"))

POISSON_RECON_BIN = os.environ.get("NPET2_POISSON_RECON_BIN", "PoissonRecon")

# API base URL for fetching profiles/landmarks when running standalone
RIBOXYZ_API_BASE = os.environ.get("NPET2_RIBOXYZ_API_URL", "http://localhost:8000")
