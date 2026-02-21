# ribctl/lib/npet2/core/run_id.py
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict


def stable_hash_dict(d: Dict[str, Any]) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _next_seq_index(runs_dir: Path) -> int:
    """
    Scan existing run directories under runs_dir for the pattern NNN_...
    and return max+1.  If none exist, returns 1.
    """
    max_idx = 0
    if runs_dir.exists():
        for d in runs_dir.iterdir():
            if d.is_dir():
                m = re.match(r"^(\d{3,})_", d.name)
                if m:
                    max_idx = max(max_idx, int(m.group(1)))
    return max_idx + 1


def compute_run_id(
    *,
    rcsb_id: str,
    pipeline_version: str,
    inputs_fp: Dict[str, str],
    config_resolved: Dict[str, Any],
    runs_dir: Path,
) -> str:
    """
    run_id = SEQ_TIMESTAMP_HASH

    Format: NNN_YYYYMMDD_HHMMSS_<hash16>
    Sequential index makes it trivial to find latest run.
    """
    blob = {
        "rcsb_id": rcsb_id.upper(),
        "pipeline_version": pipeline_version,
        "inputs": dict(sorted(inputs_fp.items())),
        "config": config_resolved,
    }
    hash_str = stable_hash_dict(blob)[:16]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    seq = _next_seq_index(runs_dir)
    return f"{seq:03d}_{timestamp}_{hash_str}"
