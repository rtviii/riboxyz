from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import shutil

from .run_id import stable_hash_dict

@dataclass(frozen=True)
class StageCacheKey:
    stage: str
    inputs_fp: dict
    params: dict
    impl_version: str = "v1"  # bump when you change semantics

    def digest(self) -> str:
        return stable_hash_dict({
            "stage": self.stage,
            "inputs_fp": self.inputs_fp,
            "params": self.params,
            "impl_version": self.impl_version,
        })[:20]

class LocalStageCache:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def entry_dir(self, key: StageCacheKey) -> Path:
        d = self.root / key.stage / key.digest()
        d.mkdir(parents=True, exist_ok=True)
        return d

    def has(self, key: StageCacheKey, required: list[str]) -> bool:
        d = self.root / key.stage / key.digest()
        if not d.exists():
            return False
        return all((d / r).exists() for r in required)

    def copy_into(self, key: StageCacheKey, dest: Path, files: list[str]) -> None:
        src = self.root / key.stage / key.digest()
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            shutil.copy2(src / f, dest / f)

    def put_from(self, key: StageCacheKey, src_dir: Path, files: list[str]) -> None:
        dst = self.entry_dir(key)
        for f in files:
            shutil.copy2(src_dir / f, dst / f)
