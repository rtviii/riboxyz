from pathlib import Path
from typing import Dict, Optional, Any

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.pipeline.setup_stage import SetupStage
from ribctl.lib.npet.pipeline.ptc_identification_stage import PTCIdentificationStage
from ribctl.lib.npet.pipeline.constriction_identification_stage import ConstrictionIdentificationStage
from ribctl.lib.npet.pipeline.exterior_mesh_stage import AlphaShapeStage
from ribctl.lib.npet.pipeline.landmark_identification_stage import LandmarkIdentificationStage
from ribctl.lib.npet.pipeline.entity_filtering_stage import EntityFilteringStage
from ribctl.lib.npet.pipeline.point_cloud_processing_stage import PointCloudProcessingStage
from ribctl.lib.npet.pipeline.clustering_stage import ClusteringStage
from ribctl.lib.npet.pipeline.refinement_stage import RefinementStage
from ribctl.lib.npet.pipeline.surface_extraction_stage import SurfaceExtractionStage
from ribctl.lib.npet.pipeline.normal_estimation_stage import NormalEstimationStage
from ribctl.lib.npet.pipeline.mesh_reconstruction_stage import MeshReconstructionStage
from ribctl.lib.npet.pipeline.validation_stage import ValidationStage
from ribctl.lib.npet.pipeline_visualization.pipeline_status_tracker import NPETProcessingTracker


def create_npet_mesh(rcsb_id: str, log_dir: Optional[Path] = None, force: bool = False) -> NPETProcessingTracker:
    """
    Creates NPET mesh for a given RCSB ID with comprehensive tracking and logging.
    
    This function orchestrates the execution of all pipeline stages, handling setup,
    error management, and tracking of the overall process.
    
    Args:
        rcsb_id: The RCSB PDB identifier
        log_dir: Directory to store logs (default: logs)
        force: Whether to force regeneration of existing assets
        
    Returns:
        NPETProcessingTracker: Processing tracker with detailed status information
    """
    print(f"Creating NPET mesh for {rcsb_id}")
    meshpath = AssetType.NPET_MESH.get_path(rcsb_id)
    if meshpath.exists():
        print(f"NPET mesh already exists for {rcsb_id}")

    rcsb_id = rcsb_id.upper()
    tracker = NPETProcessingTracker(rcsb_id, log_dir)
    
    try:
        # Initialize pipeline context to share data between stages
        context: Dict[str, Any] = {}
        
        # Define stage parameters
        alpha_shape_params = {
            "d3d_alpha": 200,
            "d3d_tol": 10,
            "d3d_offset": 3,

            "kdtree_radius": 40,
            "max_nn": 60,
            "tangent_planes_k": 20,

            "PR_depth": 6,
            "PR_ptweight": 4,
        }
        
        point_cloud_params = {
            "radius": 35,
            "height": 120,
            "voxel_size": 1,
            "atom_size": 2,
        }
        
        clustering_params = {
            "epsilon": 5.5,
            "min_samples": 600,
        }
        
        refinement_params = {
            "epsilon": 3.5,
            "min_samples": 175,
        }
        
        surface_params = {
            "alpha": 2,
            "tolerance": 1,
            "offset": 2,
        }
        
        normal_params = {
            "kdtree_radius": 10,
            "kdtree_max_nn": 15,
            "correction_tangent_planes_n": 10,
        }
        
        mesh_params = {
            "depth": 6,
            "ptweight": 3,
        }
        
        # Run the Setup stage
        setup_stage = SetupStage(rcsb_id, tracker, None, force)
        setup_result = setup_stage.run(context)
        context.update(setup_result)
        
        # Create artifacts directory
        artifacts_dir = Path(context["ro"].assets.paths.dir) / "artifacts"
        artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        # Run the PTC identification stage
        ptc_stage = PTCIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        ptc_result = ptc_stage.run(context)
        context.update(ptc_result)
        
        # Run the constriction identification stage
        constriction_stage = ConstrictionIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        constriction_result = constriction_stage.run(context)
        context.update(constriction_result)
        
        # Run the alpha shape stage
        alpha_stage = AlphaShapeStage(rcsb_id, tracker, artifacts_dir, True, alpha_shape_params)
        alpha_result = alpha_stage.run(context)
        context.update(alpha_result)
        
        # Run the landmark identification stage (verifies PTC and constriction)
        landmark_stage = LandmarkIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        landmark_result = landmark_stage.run(context)
        context.update(landmark_result)
        
        # Run the entity filtering stage
        entity_filtering_stage = EntityFilteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"]
        )
        entity_result = entity_filtering_stage.run(context)
        context.update(entity_result)
        
        # Run the point cloud processing stage
        point_cloud_stage = PointCloudProcessingStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"],
            voxel_size=point_cloud_params["voxel_size"],
            atom_size=point_cloud_params["atom_size"]
        )
        point_cloud_result = point_cloud_stage.run(context)
        context.update(point_cloud_result)
        
        # Run the clustering stage
        clustering_stage = ClusteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=clustering_params["epsilon"],
            min_samples=clustering_params["min_samples"]
        )
        clustering_result = clustering_stage.run(context)
        context.update(clustering_result)
        
        # Run the refinement stage
        refinement_stage = RefinementStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=refinement_params["epsilon"],
            min_samples=refinement_params["min_samples"]
        )
        refinement_result = refinement_stage.run(context)
        context.update(refinement_result)
        
        # Run the surface extraction stage
        surface_stage = SurfaceExtractionStage(
            rcsb_id, tracker, artifacts_dir, force,
            alpha=surface_params["alpha"],
            tolerance=surface_params["tolerance"],
            offset=surface_params["offset"]
        )
        surface_result = surface_stage.run(context)
        context.update(surface_result)
        
        # Run the normal estimation stage
        normal_stage = NormalEstimationStage(
            rcsb_id, tracker, artifacts_dir, force,
            kdtree_radius=normal_params["kdtree_radius"],
            kdtree_max_nn=normal_params["kdtree_max_nn"],
            correction_tangent_planes_n=normal_params["correction_tangent_planes_n"]
        )
        normal_result = normal_stage.run(context)
        context.update(normal_result)
        
        # Run the mesh reconstruction stage
        mesh_stage = MeshReconstructionStage(
            rcsb_id, tracker, artifacts_dir, force,
            depth=mesh_params["depth"],
            ptweight=mesh_params["ptweight"]
        )
        mesh_result = mesh_stage.run(context)
        context.update(mesh_result)
        
        # Run the validation stage
        validation_stage = ValidationStage(rcsb_id, tracker, artifacts_dir, force)
        validation_result = validation_stage.run(context)
        context.update(validation_result)
        
        # Complete pipeline and set final status
        tracker.complete_processing(True, watertight=context.get("watertight", False))
        
    except Exception as e:
        # Handle any uncaught exceptions
        if tracker.current_stage:
            tracker.end_stage(tracker.current_stage, False, error=e)
        tracker.complete_processing(False)
    
    return tracker