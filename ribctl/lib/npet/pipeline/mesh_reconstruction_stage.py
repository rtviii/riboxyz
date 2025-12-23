from pathlib import Path
from typing import Any, Dict
import os

from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage


class MeshReconstructionStage(NPETPipelineStage):
    """
    Stage for reconstructing the mesh from the point cloud with normals.
    
    This stage applies Poisson surface reconstruction to create a 
    3D mesh representing the tunnel surface.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 depth: int = 6,
                 ptweight: int = 3):
        """
        Initialize the mesh reconstruction stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            depth: Depth parameter for Poisson reconstruction
            ptweight: Point weight for Poisson reconstruction
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.depth = depth
        self.ptweight = ptweight
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.MESH_RECONSTRUCTION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "depth": self.depth,
            "ptweight": self.ptweight,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Reconstruct the mesh from the point cloud with normals.
        
        Steps:
        1. Apply Poisson reconstruction to create the mesh
        2. Save the mesh and related files
        """
        normals_pcd_path = context["normals_pcd_path"]
        meshpath = context["meshpath"]
        
        # Apply Poisson reconstruction
        apply_poisson_reconstruction(
            str(normals_pcd_path),
            meshpath,
            recon_depth=self.depth,
            recon_pt_weight=self.ptweight,
        )
        
        # Add mesh as artifact
        self.tracker.add_artifact(self.stage, Path(meshpath))
        
        # Also save ASCII version as artifact if it exists
        ascii_path = Path(str(meshpath).split(".")[0] + "_ascii.ply")
        if ascii_path.exists():
            self.tracker.add_artifact(self.stage, ascii_path)
        
        return {
            "meshpath": meshpath
        }