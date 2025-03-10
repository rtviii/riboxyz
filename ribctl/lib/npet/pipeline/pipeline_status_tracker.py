import json
import time
import traceback
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List

class ProcessingStage(str, Enum):
    """Enum defining the stages of the NPET mesh pipeline"""
    SETUP = "setup"
    LANDMARK_IDENTIFICATION = "landmark_identification"
    ENTITY_FILTERING = "entity_filtering"
    POINT_CLOUD_PROCESSING = "point_cloud_processing"
    CLUSTERING = "clustering"
    REFINEMENT = "refinement"
    SURFACE_EXTRACTION = "surface_extraction"
    NORMAL_ESTIMATION = "normal_estimation"
    MESH_RECONSTRUCTION = "mesh_reconstruction"
    VALIDATION = "validation"
    COMPLETE = "complete"

class ProcessingStatus(str, Enum):
    """Status of a processing stage"""
    PENDING = "pending"
    IN_PROGRESS = "in_progress"
    SUCCESS = "success"
    FAILURE = "failure"
    SKIPPED = "skipped"

class NPETProcessingTracker:
    """Tracks the processing status of NPET mesh generation for a structure"""
    
    def __init__(self, rcsb_id: str, output_dir: Optional[Path] = None):
        self.rcsb_id = rcsb_id.upper()
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.current_stage: ProcessingStage = ProcessingStage.SETUP
        self.output_dir = output_dir or Path(f"logs")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize stage tracking
        self.stages: Dict[ProcessingStage, Dict[str, Any]] = {
            stage: {
                "status": ProcessingStatus.PENDING,
                "start_time": None,
                "end_time": None,
                "duration": None,
                "error": None,
                "artifacts": []
            } for stage in ProcessingStage
        }
        
        # Initialize overall status
        self.status: Dict[str, Any] = {
            "rcsb_id": self.rcsb_id,
            "overall_status": ProcessingStatus.PENDING,
            "start_time": self.start_time,
            "end_time": None,
            "duration": None,
            "stages": self.stages,
            "summary": {
                "success": False,
                "failed_stage": None,
                "error_summary": None,
                "watertight": False,
                "artifacts_generated": []
            }
        }
        
        # Save initial status
        self._save_status()
    
    def begin_stage(self, stage: ProcessingStage) -> None:
        """Mark the beginning of a processing stage"""
        self.current_stage = stage
        now = time.time()
        self.stages[stage].update({
            "status": ProcessingStatus.IN_PROGRESS,
            "start_time": now
        })
        self._save_status()
    
    def end_stage(self, stage: ProcessingStage, success: bool, 
                 artifacts: Optional[List[Path]] = None,
                 error: Optional[Exception] = None) -> None:
        """Mark the end of a processing stage"""
        now = time.time()
        start_time = self.stages[stage]["start_time"] or now
        
        status = ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE
        
        error_info = None
        if error:
            error_info = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc()
            }
        
        self.stages[stage].update({
            "status": status,
            "end_time": now,
            "duration": now - start_time,
            "error": error_info,
            "artifacts": [str(a) for a in (artifacts or [])]
        })
        
        if not success:
            # Update summary with failure info
            self.status["summary"]["success"] = False
            self.status["summary"]["failed_stage"] = stage
            self.status["summary"]["error_summary"] = str(error) if error else "Unknown error"
            
            # Mark remaining stages as skipped
            all_stages = list(ProcessingStage)
            current_index = all_stages.index(stage)
            for skip_stage in all_stages[current_index+1:]:
                if skip_stage != ProcessingStage.COMPLETE:
                    self.stages[skip_stage]["status"] = ProcessingStatus.SKIPPED
        
        self._save_status()
    
    def complete_processing(self, success: bool, watertight: bool = False) -> None:
        """Mark the overall processing as complete"""
        now = time.time()
        self.end_time = now
        
        # Get all generated artifacts
        artifacts = []
        for stage_info in self.stages.values():
            artifacts.extend(stage_info.get("artifacts", []))
        
        self.status.update({
            "overall_status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "end_time": now,
            "duration": now - self.start_time,
            "summary": {
                "success": success,
                "watertight": watertight,
                "artifacts_generated": artifacts
            }
        })
        
        # Update the complete stage
        self.stages[ProcessingStage.COMPLETE].update({
            "status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "start_time": now,
            "end_time": now,
            "duration": 0
        })
        
        self._save_status()
    
    def add_artifact(self, stage: ProcessingStage, artifact_path: Path) -> None:
        """Add an artifact to the specified stage"""
        if str(artifact_path) not in self.stages[stage].get("artifacts", []):
            self.stages[stage].setdefault("artifacts", []).append(str(artifact_path))
            self._save_status()
    
    def _save_status(self) -> None:
        """Save the current status to a JSON file"""
        log_file = self.output_dir / f"{self.rcsb_id}_processing_log.json"
        with open(log_file, 'w') as f:
            json.dump(self.status, f, indent=2, default=str)