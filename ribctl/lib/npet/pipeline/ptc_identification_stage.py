import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import ProcessingStage


class PTCIdentificationStage(NPETPipelineStage):
    """
    Stage for identifying the PTC (Peptidyl Transferase Center) location.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.PTC_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # No configurable parameters for this stage
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Identify the PTC location based on the ribosome structure.
        """
        try:
            # Get PTC location
            ptc_info = PTC_location(self.rcsb_id)
            ptc_pt = np.array(ptc_info.location)
            
            # Save PTC info as JSON
            ptc_json_path = AssetType.PTC.get_path(self.rcsb_id)
            with open(ptc_json_path, 'w') as f:
                f.write(ptc_info.model_dump_json())
            
            # Track the artifact
            self.tracker.add_artifact(self.stage, ptc_json_path)
            
            return {
                "ptc_info": ptc_info,
                "ptc_pt": ptc_pt
            }
        except Exception as e:
            raise RuntimeError(f"Failed to identify PTC: {str(e)}") from e