import json
import numpy as np
from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import ProcessingStage
from ribctl.lib.schema.types_ribosome import ConstrictionSite


class ConstrictionIdentificationStage(NPETPipelineStage):
    """
    Stage for identifying the constriction site in the ribosome tunnel.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.CONSTRICTION_IDENTIFICATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # No configurable parameters for this stage
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Identify the constriction site within the ribosome tunnel.
        """
        try:
            # Get constriction point
            constriction_pt = get_constriction(self.rcsb_id)
            
            # Save constriction point as JSON
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(self.rcsb_id)
            with open(constriction_path, 'w') as f:
                json.dump(ConstrictionSite(location=constriction_pt.tolist()).model_dump(), f)
            
            # Track the artifact
            self.tracker.add_artifact(self.stage, constriction_path)
            
            return {
                "constriction_pt": constriction_pt
            }
        except Exception as e:
            raise RuntimeError(f"Failed to identify constriction site: {str(e)}") from e