import os
from pathlib import Path
from typing import Any, Dict
from ribctl.lib.npet.alphalib import validate_mesh_pyvista
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_status_tracker import ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_mesh


class ValidationStage(NPETPipelineStage):
    """
    Stage for validating the final mesh.
    
    This stage checks if the generated mesh is watertight and meets
    quality standards. It also generates visualizations of the final mesh.
    """
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.VALIDATION
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {}
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Validate the final mesh and ensure it meets quality standards.
        
        Steps:
        1. Generate visualization of the final mesh
        2. Check if the mesh is watertight
        3. If not watertight, remove the mesh files and raise an error
        """
        meshpath = context["meshpath"]
        
        # Save mesh visualization if possible
        try:
            mesh_viz_path = self.artifacts_dir / f"{self.rcsb_id}_mesh.png"
            visualize_mesh(meshpath, self.rcsb_id, output_path=str(mesh_viz_path))
            self.tracker.add_artifact(self.stage, mesh_viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save mesh visualization: {str(viz_error)}")
        
        # Check watertightness
        watertight = validate_mesh_pyvista(meshpath)
        
        if not watertight:
            error_msg = "Mesh is not watertight"
            print("XXXX Watertightness check failed, removing", meshpath, " XXXX")
            os.remove(meshpath)
            
            # Also remove ASCII version if it exists
            ascii_path = str(meshpath).split(".")[0] + "_ascii.ply"
            if os.path.exists(ascii_path):
                os.remove(ascii_path)
                
            raise ValueError(error_msg)
        
        return {
            "watertight": watertight
        }