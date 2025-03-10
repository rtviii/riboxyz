import copy
import json
import os
import time
import traceback
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List

class ProcessingStage(str, Enum):
    """Enum defining the stages of the NPET mesh pipeline"""
    SETUP = "setup"
    PTC_IDENTIFICATION = "ptc_identification"           # New stage
    CONSTRICTION_IDENTIFICATION = "constriction_identification"  # New stage
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
        
        # Print stage start
        print(f"\n⏳ Starting stage: {stage}")
        
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
                    print(f"⏩ Stage {stage} skipped (parameters unchanged, artifacts exist)")
                    
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
                    
            # Print failure message with emoji
            print(f"❌ Stage FAILED: {stage} - {str(error) if error else 'Unknown error'}")
        else:
            # Print success message with emoji
            duration = now - start_time
            print(f"✅ Stage completed: {stage} ({duration:.2f}s)")
        
        # Print artifact count
        artifact_count = len(self.stages[stage].get("artifacts", []))
        if artifact_count > 0:
            print(f"   {artifact_count} artifacts generated")
        
        self._save_status()
    
    def complete_processing(self, success: bool, watertight: bool = False) -> None:
        """Mark the overall processing as complete"""
        now = time.time()
        self.end_time = now
        total_duration = now - self.start_time
        
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
            "duration": total_duration,
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
        
        # Don't print summary here - let the caller handle it
        # self.print_summary()
        
        self._save_status()
    
    def add_artifact(self, stage: ProcessingStage, artifact_path: Path) -> None:
        """Add an artifact to the specified stage"""
        artifact_str = str(artifact_path)
        
        if artifact_str not in self.stages[stage].get("artifacts", []):
            self.stages[stage].setdefault("artifacts", []).append(artifact_str)
            self._save_status()
    
    def print_summary(self) -> None:
        """Print a summary of all stages with their status"""
        print("\n" + "=" * 60)
        print(f"PROCESSING SUMMARY FOR {self.rcsb_id}")
        print("=" * 60)
        
        # Calculate the widest stage name for nice formatting
        stage_width = max(len(str(stage)) for stage in ProcessingStage) + 2
        
        # Print status for each stage
        for stage in ProcessingStage:
            if stage != ProcessingStage.COMPLETE:
                status = self.stages[stage]["status"]
                
                # Select appropriate emoji and color based on status
                emoji = {
                    ProcessingStatus.SUCCESS: "✅",
                    ProcessingStatus.FAILURE: "❌",
                    ProcessingStatus.SKIPPED: "⏩",
                    ProcessingStatus.IN_PROGRESS: "⏳",
                    ProcessingStatus.PENDING: "⏱️",
                }.get(status, "❓")
                
                # Get duration if available
                duration_str = ""
                if self.stages[stage]["duration"] is not None:
                    duration = self.stages[stage]["duration"]
                    if duration < 0.001:  # Skipped stages
                        duration_str = "(skipped)"
                    else:
                        duration_str = f"({duration:.2f}s)"
                
                # Get artifacts if available
                artifact_count = len(self.stages[stage].get("artifacts", []))
                artifact_str = f"{artifact_count} artifacts" if artifact_count > 0 else ""
                
                # Get error if available
                error_str = ""
                if self.stages[stage]["error"]:
                    error_str = f"ERROR: {self.stages[stage]['error']['message']}"
                
                # Format the line
                stage_name = f"{stage}".ljust(stage_width)
                status_name = f"{status.name}".ljust(10)
                duration_str = duration_str.ljust(15)
                artifact_str = artifact_str.ljust(15)
                
                print(f"{emoji} {stage_name} {status_name} {duration_str} {artifact_str} {error_str}")
        
        print("-" * 60)
        
        # Print overall status
        total_duration = self.status.get("duration", 0)
        duration_str = f"Total time: {total_duration:.2f}s" if total_duration else ""
        
        success = self.status["summary"]["success"]
        if success:
            print(f"✅ PIPELINE COMPLETED SUCCESSFULLY  {duration_str}")
            if self.status["summary"]["watertight"]:
                mesh_files = [a for a in self.status["summary"]["artifacts_generated"] 
                            if a.endswith('.ply') and not a.endswith('_normal_estimated_pcd.ply')]
                if mesh_files:
                    print(f"   Watertight mesh: {os.path.basename(mesh_files[0])}")
        else:
            failed_stage = self.status["summary"]["failed_stage"]
            error = self.status["summary"]["error_summary"]
            print(f"❌ PIPELINE FAILED at stage {failed_stage}  {duration_str}")
            if error:
                print(f"   Error: {error}")
        
        print("=" * 60)
    
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