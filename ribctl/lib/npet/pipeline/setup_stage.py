from pathlib import Path
from typing import Any, Dict

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import ProcessingStage
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager
from ribctl.ribosome_ops import RibosomeOps


class SetupStage(NPETPipelineStage):
    """
    Setup stage - initializes resources and validates inputs.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.SETUP
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        # This stage doesn't have configurable parameters
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Initialize the pipeline by setting up paths and resources.
        """
        # Initialize asset managers
        assets = TunnelMeshAssetsManager(self.rcsb_id)
        ro = RibosomeOps(self.rcsb_id)
        
        # Set up paths
        cifpath = AssetType.MMCIF.get_path(self.rcsb_id)
        ashapepath = AssetType.ALPHA_SHAPE.get_path(self.rcsb_id)
        meshpath = AssetType.NPET_MESH.get_path(self.rcsb_id)
        
        # Track the CIF file as an artifact
        self.tracker.add_artifact(self.stage, Path(cifpath))
        
        return {
            "assets": assets,
            "ro": ro,
            "profile": ro.profile,
            "cifpath": cifpath,
            "ashapepath": ashapepath,
            "meshpath": meshpath
        }