# npet2/run.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

from npet2.adapters.standalone_providers import FileStructureProvider, FileLandmarkProvider
from npet2.core.config import RunConfig, SETTINGS
from npet2.core.manifest import RunManifest
from npet2.core.run_id import compute_run_id
from npet2.core.store import LocalRunStore
from npet2.core.types import StageContext
from npet2.core.pipeline import Pipeline
from npet2.core.cache import LocalStageCache

from npet2.stages.bootstrap import Stage00Inputs, Stage10Landmarks
from npet2.stages.grid_refine import Stage55GridRefine
from npet2.stages.legacy_minimal import (
    Stage20ExteriorShell,
    Stage30RegionAtoms,
    Stage40EmptySpace,
    Stage50Clustering,
    Stage70MeshValidate,
)


def _pipeline_version() -> str:
    return "npet2-dev"


def run_npet2(
    rcsb_id: str,
    config: Optional[RunConfig] = None,
    *,
    structure_provider=None,
    landmark_provider=None,
) -> StageContext:
    rcsb_id = rcsb_id.upper()
    config = config or RunConfig()

    if structure_provider is None or landmark_provider is None:
        raise ValueError(
            "Both structure_provider and landmark_provider are required. "
            "Use FileStructureProvider / FileLandmarkProvider for standalone mode."
        )

    config_resolved = asdict(config)
    inputs_fp = {
        "structure": structure_provider.fingerprint(rcsb_id),
        "landmarks": landmark_provider.fingerprint(rcsb_id),
    }

    struct_runs_dir = SETTINGS.runs_root / rcsb_id
    struct_runs_dir.mkdir(parents=True, exist_ok=True)

    run_id = compute_run_id(
        rcsb_id=rcsb_id,
        pipeline_version=_pipeline_version(),
        inputs_fp=inputs_fp,
        config_resolved=config_resolved,
        runs_dir=struct_runs_dir,
    )

    run_dir = struct_runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    manifest = RunManifest(
        rcsb_id=rcsb_id,
        run_id=run_id,
        pipeline_version=_pipeline_version(),
        inputs={"fingerprints": inputs_fp},
        config_resolved=config_resolved,
    )
    store = LocalRunStore(run_dir=run_dir, manifest=manifest)

    ctx = StageContext(
        run_id=run_id,
        rcsb_id=rcsb_id,
        config=config,
        store=store,
        inputs={
            "structure_provider": structure_provider,
            "landmark_provider": landmark_provider,
            "stage_cache": LocalStageCache(SETTINGS.cache_root),
            "inputs_fp": inputs_fp,
        },
    )

    pipeline = Pipeline([
        Stage00Inputs(),
        Stage10Landmarks(),
        Stage20ExteriorShell(),
        Stage30RegionAtoms(),
        Stage40EmptySpace(),
        Stage50Clustering(),
        Stage55GridRefine(),
        Stage70MeshValidate(),
    ])

    pipeline.run(ctx)
    return ctx
