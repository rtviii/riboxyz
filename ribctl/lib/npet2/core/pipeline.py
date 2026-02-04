# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict
from typing import Any, Dict, List, Optional

from .types import StageContext


class Stage(ABC):
    """
    A stage consumes/produces via ctx.inputs and ctx.artifacts.
    Stage key is a stable string like "10_landmarks".
    """

    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]:
        ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None:
        ...


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        for stage in self.stages:
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)
            try:
                stage.run(ctx)
                ctx.store.end_stage(stage.key, success=True)
            except Exception as e:
                ctx.store.end_stage(stage.key, success=False, note=str(e))
                ctx.store.finalize(success=False, error=str(e))
                raise
        ctx.store.finalize(success=True)
        return ctx
