# ribctl/lib/npet2/core/run_id.py
from __future__ import annotations

import hashlib
import json
from datetime import datetime
from typing import Any, Dict


def stable_hash_dict(d: Dict[str, Any]) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def compute_run_id(
    *, 
    rcsb_id: str, 
    pipeline_version: str, 
    inputs_fp: Dict[str, str], 
    config_resolved: Dict[str, Any]
) -> str:
    """
    run_id = timestamp_hash
    
    Format: YYYYMMDD_HHMMSS_<hash16>
    This allows chronological sorting while keeping collision resistance.
    """
    blob = {
        "rcsb_id": rcsb_id.upper(),
        "pipeline_version": pipeline_version,
        "inputs": dict(sorted(inputs_fp.items())),
        "config": config_resolved,
    }
    hash_str = stable_hash_dict(blob)[:16]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{timestamp}_{hash_str}"
