from abc import ABC, abstractmethod
from pathlib import Path
from typing import Any, Dict, Optional

from ribctl.lib.npet.visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage



class NPETPipelineStage(ABC):
    """
    Base class for all NPET pipeline stages.
    
    Each stage handles a specific part of the mesh creation process,
    with standardized interfaces for input/output and error handling.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False):
        """
        Initialize a pipeline stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker for logging progress
            artifacts_dir: Directory to store output artifacts
            force: Whether to force regeneration even if artifacts exist
        """
        self.rcsb_id = rcsb_id.upper()
        self.tracker = tracker
        self.artifacts_dir = artifacts_dir
        self.force = force
    
    @property
    @abstractmethod
    def stage(self) -> ProcessingStage:
        """The processing stage this class handles."""
        pass
    
    @property
    @abstractmethod
    def stage_params(self) -> Dict[str, Any]:
        """Parameters used by this stage for tracking and reproducibility."""
        pass
    
    @abstractmethod
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Execute the processing for this stage.
        
        Args:
            context: Input data from previous stages
            
        Returns:
            Dict with outputs for next stages
        """
        pass
    
    def run(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run the stage with progress tracking and error handling.
        
        Args:
            context: Input data from previous stages
            
        Returns:
            Dict with outputs for next stages
        """
        should_process = self.tracker.begin_stage(self.stage, self.stage_params)
        
        if not should_process and not self.force:
            print(f"Skipping {self.stage} (artifacts exist and parameters unchanged)")
            return {}
            
        try:
            results = self.process(context)
            self.tracker.end_stage(self.stage, True)
            return results
        except Exception as e:
            self.tracker.end_stage(self.stage, False, error=e)
            # Re-raise to stop pipeline
            raise