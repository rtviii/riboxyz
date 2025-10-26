import json
import numpy as np
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.visualization.pipeline_status_tracker import ProcessingStage
from ribctl.lib.schema.types_ribosome import PTCInfo, ConstrictionSite


class LandmarkIdentificationStage(NPETPipelineStage):
    """
    Stage for loading and verifying landmark points (PTC and constriction site).
    
    This stage serves as a verification point to ensure both the PTC and
    constriction site are available for subsequent stages.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.LANDMARK_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # This stage doesn't have configurable parameters
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Load and verify landmark points from previous stages or files.
        
        The stage either uses the PTC and constriction points from previous
        stages or loads them from disk if they weren't already computed.
        """
        try:
            # Check if we already have the landmarks from previous stages
            ptc_pt = context.get("ptc_pt")
            constriction_pt = context.get("constriction_pt")
            
            # If not available in context, try to load from files
            if ptc_pt is None or constriction_pt is None:
                # Load PTC and constriction points from previous stages
                ptc_path = AssetType.PTC.get_path(self.rcsb_id)
                constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
                
                # If they exist, use them directly
                if ptc_path.exists() and constriction_path.exists():
                    with open(ptc_path, 'r') as f:
                        ptc_info = PTCInfo.model_validate(json.load(f))
                    with open(constriction_path, 'r') as f:
                        constriction_site = ConstrictionSite.model_validate(json.load(f))
                    ptc_pt = np.array(ptc_info.location)
                    constriction_pt = np.array(constriction_site.location)
                else:
                    raise FileNotFoundError("PTC or constriction site not found")
            
            # Add existing files as artifacts if they weren't already tracked
            ptc_path = AssetType.PTC.get_path(self.rcsb_id)
            if ptc_path.exists():
                self.tracker.add_artifact(self.stage, ptc_path)
                
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
            if constriction_path.exists():
                self.tracker.add_artifact(self.stage, constriction_path)
    
            return {
                "ptc_pt": ptc_pt,
                "constriction_pt": constriction_pt
            }
            
        except Exception as e:
            raise RuntimeError(f"Failed to identify landmarks: {str(e)}") from e