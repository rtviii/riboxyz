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


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        for stage in self.stages:
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)

            t0 = time.perf_counter()
            print(f"[npet2] >>> {stage.key} start")

            try:
                stage.run(ctx)
                dt = time.perf_counter() - t0
                print(f"[npet2] <<< {stage.key} done in {dt:,.2f}s")
                ctx.store.end_stage(stage.key, success=True, note=f"elapsed_s={dt:.3f}")
            except Exception as e:
                dt = time.perf_counter() - t0
                print(f"[npet2] !!! {stage.key} FAILED after {dt:,.2f}s: {e}")
                ctx.store.end_stage(
                    stage.key, success=False, note=f"elapsed_s={dt:.3f} err={e}"
                )
                ctx.store.finalize(success=False, error=str(e))
                raise

        ctx.store.finalize(success=True)
        return ctx
