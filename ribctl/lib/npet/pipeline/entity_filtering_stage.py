import numpy as np
from pathlib import Path
from typing import Any, Dict, List

from ribctl.lib.npet.kdtree_approach import (
    ribosome_entities,
    filter_residues_parallel,
)
from ribctl.lib.npet.pipeline.base_stage import NPETPipelineStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.various_visualization import visualize_filtered_residues


class EntityFilteringStage(NPETPipelineStage):
    """
    Stage for filtering ribosome entities to only those within the tunnel region.
    """
    
    def __init__(self, 
                 rcsb_id: str, 
                 tracker: NPETProcessingTracker,
                 artifacts_dir: Path,
                 force: bool = False,
                 radius: float = 35,
                 height: float = 120):
        """
        Initialize the entity filtering stage.
        
        Args:
            rcsb_id: The RCSB PDB identifier
            tracker: Processing tracker
            artifacts_dir: Directory for artifacts
            force: Whether to force regeneration
            radius: Radius of the tunnel cylinder
            height: Height of the tunnel cylinder
        """
        super().__init__(rcsb_id, tracker, artifacts_dir, force)
        self.radius = radius
        self.height = height
    
    @property
    def stage(self) -> ProcessingStage:
        return ProcessingStage.ENTITY_FILTERING
    
    @property
    def stage_params(self) -> Dict[str, Any]:
        return {
            "radius": self.radius,
            "height": self.height,
        }
    
    def process(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """
        Filter ribosome entities to only those relevant to the tunnel.
        """
        cifpath = context["cifpath"]
        ptc_pt = context["ptc_pt"]
        constriction_pt = context["constriction_pt"]
        profile = context["profile"]
        
        # Determine tunnel debris chains to exclude
        tunnel_debris = {
            "3J7Z": ["a", "7"],
            "5GAK": ["z"],
            "5NWY": ["s"],
            "7A5G": ["Y2"],
        }
        
        # For mitochondrial ribosomes, also exclude mL45 chain
        if profile.mitochondrial:
            try:
                chain = context["ro"].get_poly_by_polyclass("mL45")
                tunnel_debris[self.rcsb_id] = [chain.auth_asym_id]
            except Exception as chain_error:
                print("Mitochondrial mL45 chain not found:", str(chain_error))
        
        # Get all ribosome entities
        residues = ribosome_entities(
            self.rcsb_id,
            cifpath,
            "R",
            tunnel_debris[self.rcsb_id] if self.rcsb_id in tunnel_debris else [],
        )
        
        # Filter to only those within tunnel region
        filtered_residues = filter_residues_parallel(
            residues, ptc_pt, constriction_pt, self.radius, self.height
        )
        
        # Extract atom coordinates
        filtered_points = np.array([
            atom.get_coord()
            for residue in filtered_residues
            for atom in residue.child_list
        ])
        
        # Save filtered points
        filtered_points_path = self.artifacts_dir / f"{self.rcsb_id}_filtered_points.npy"
        np.save(filtered_points_path, filtered_points)
        self.tracker.add_artifact(self.stage, filtered_points_path)
        
        # Save visualization if possible
        try:
            viz_path = self.artifacts_dir / f"{self.rcsb_id}_filtered_residues.png"
            visualize_filtered_residues(
                filtered_residues, residues, ptc_pt, constriction_pt, self.radius, self.height,
                output_path=str(viz_path)
            )
            self.tracker.add_artifact(self.stage, viz_path)
        except Exception as viz_error:
            print(f"Warning: Could not save visualization: {str(viz_error)}")
        
        return {
            "filtered_residues": filtered_residues,
            "filtered_points": filtered_points
        }