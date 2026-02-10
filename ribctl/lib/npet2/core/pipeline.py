# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List
import time

from .types import StageContext


class Stage(ABC):
    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]: ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None: ...


def _fmt_value(v: Any) -> str:
    """Compact display for a parameter value."""
    if isinstance(v, float):
        # drop trailing zeros but keep one decimal
        return f"{v:g}"
    if isinstance(v, list) and len(v) > 3:
        return f"[{len(v)} items]"
    return str(v)


def _log_params(params: Dict[str, Any], prefix: str) -> None:
    if not params:
        return
    items = [f"{k}={_fmt_value(v)}" for k, v in params.items()
             if not isinstance(v, (dict, list))]
    # nested dicts/lists get their own lines
    nested = {k: v for k, v in params.items() if isinstance(v, (dict, list))}

    if items:
        line = ", ".join(items)
        # wrap at ~100 chars
        if len(line) > 100:
            mid = len(items) // 2
            print(f"  [{prefix}] {', '.join(items[:mid])}")
            print(f"  [{prefix}] {', '.join(items[mid:])}")
        else:
            print(f"  [{prefix}] {line}")
    for k, v in nested.items():
        if isinstance(v, list) and all(isinstance(x, dict) for x in v):
            for i, entry in enumerate(v):
                sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in entry.items())
                print(f"  [{prefix}] {k}[{i}]: {sub}")
        elif isinstance(v, dict):
            sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in v.items())
            print(f"  [{prefix}] {k}: {sub}")


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        n = len(self.stages)
        wall = 60

        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | run {ctx.run_id}")
        print(f"  stages: {n} | config: {type(ctx.config).__name__}")
        print("=" * wall)

        t_total = time.perf_counter()

        for i, stage in enumerate(self.stages, 1):
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)

            print()
            print(f"--- [{i}/{n}] {stage.key} " + "-" * max(0, wall - len(stage.key) - 12))
            _log_params(params, stage.key)

            t0 = time.perf_counter()
            try:
                stage.run(ctx)
                dt = time.perf_counter() - t0
                print(f"  [{stage.key}] done in {dt:,.2f}s")
                ctx.store.end_stage(stage.key, success=True, note=f"elapsed_s={dt:.3f}")
            except Exception as e:
                dt = time.perf_counter() - t0
                print(f"  [{stage.key}] FAILED after {dt:,.2f}s: {e}")
                ctx.store.end_stage(
                    stage.key, success=False, note=f"elapsed_s={dt:.3f} err={e}"
                )
                ctx.store.finalize(success=False, error=str(e))
                raise

        dt_total = time.perf_counter() - t_total
        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | completed in {dt_total:,.2f}s")
        print(f"  run_dir: {ctx.store.run_dir}")
        print("=" * wall)
        print()

        ctx.store.finalize(success=True)
        return ctx