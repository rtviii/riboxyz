import copy
import json
import time
import traceback
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List

class ProcessingStage(str, Enum):
    """Enum defining the stages of the NPET mesh pipeline"""
    SETUP = "setup"
    ALPHA_SHAPE = "alpha_shape"
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
                "artifacts": [],
                "parameters": {}  # New field to store parameters
            } for stage in ProcessingStage
        }
        
        # Initialize overall status
        self.status: Dict[str, Any] = {
            "rcsb_id": self.rcsb_id,
            "overall_status": ProcessingStatus.PENDING,
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
    
    def begin_stage(self, stage: ProcessingStage, parameters: Optional[Dict[str, Any]] = None) -> bool:
        """
        Mark the beginning of a processing stage.
        
        Args:
            stage: The processing stage
            parameters: Optional dictionary of parameters for this stage
            
        Returns:
            bool: True if the stage should be processed, False if it can be skipped
        """
        self.current_stage = stage
        now = time.time()
        
        # Check if we can skip processing based on parameters
        can_skip = False
        if parameters is not None:
            # Get previous parameters and artifacts
            prev_params = self.stages[stage].get("parameters", {})
            prev_artifacts = self.stages[stage].get("artifacts", [])
            
            # Check if parameters match and artifacts exist
            if prev_params and prev_artifacts and prev_params == parameters:
                # Check that artifacts actually exist on disk
                artifacts_exist = all(os.path.exists(artifact) for artifact in prev_artifacts)
                if artifacts_exist:
                    can_skip = True
                    print(f"Stage {stage} can be skipped (parameters unchanged, artifacts exist)")
                    
                    # Update status but keep existing data
                    self.stages[stage].update({
                        "status": ProcessingStatus.SUCCESS,
                        "start_time": now,
                        "end_time": now,
                        "duration": 0.0,  # Zero duration for skipped stages
                    })
                    self._save_status()
                    return False  # Skip processing
        
        # Proceed with processing
        self.stages[stage].update({
            "status": ProcessingStatus.IN_PROGRESS,
            "start_time": now,
            "parameters": parameters or {}  # Store new parameters
        })
        self._save_status()
        return True  # Process the stage
    
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
        
        # Create a dict with the updates, but don't include artifacts unless provided
        update_dict = {
            "status": status,
            "end_time": now,
            "duration": now - start_time,
            "error": error_info,
        }
        
        # Only update the artifacts list if explicitly provided
        if artifacts is not None:
            update_dict["artifacts"] = [str(a) for a in artifacts]
        
        # Apply the updates
        self.stages[stage].update(update_dict)
        
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
        
        # Print artifact count for debugging
        print(f"Stage {stage} now has {len(self.stages[stage].get('artifacts', []))} artifacts")
        
        self._save_status()
    
    def complete_processing(self, success: bool, watertight: bool = False) -> None:
        """Mark the overall processing as complete"""
        now = time.time()
        self.end_time = now
        
        # Get all generated artifacts with debugging
        artifacts = []
        for stage_name, stage_info in self.stages.items():
            stage_artifacts = stage_info.get("artifacts", [])
            if stage_artifacts:
                print(f"Found {len(stage_artifacts)} artifacts in stage {stage_name}: {stage_artifacts}")
                artifacts.extend(stage_artifacts)
        
        print(f"Total artifacts collected: {len(artifacts)}")
        
        self.status.update({
            "overall_status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
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
            "duration": 0
        })
        
        self._save_status()
    
    def add_artifact(self, stage: ProcessingStage, artifact_path: Path) -> None:
        """Add an artifact to the specified stage"""
        artifact_str = str(artifact_path)
        print(f"Adding artifact to stage {stage}: {artifact_str}")
        print(f"Artifact exists: {artifact_path.exists()}")
        
        if artifact_str not in self.stages[stage].get("artifacts", []):
            print(f"  - Artifact not already in list, adding it")
            self.stages[stage].setdefault("artifacts", []).append(artifact_str)
            
            # Verify after adding
            print(f"  - Current artifacts in stage: {self.stages[stage]['artifacts']}")
            self._save_status()
        else:
            print(f"  - Artifact already in list, skipping")
    
    def _save_status(self) -> None:
        """Save the current status to a JSON file"""
        log_file = self.output_dir / f"{self.rcsb_id}_processing_log.json"
        
        # Create a simplified version of the status for JSON output
        output_status = copy.deepcopy(self.status)
        
        # Simplify the time representation - remove start_time and end_time
        for stage_name, stage_info in output_status["stages"].items():
            if "start_time" in stage_info:
                del stage_info["start_time"]
            if "end_time" in stage_info:
                del stage_info["end_time"]
        
        with open(log_file, 'w') as f:
            json.dump(output_status, f, indent=2, default=str)