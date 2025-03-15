import json
import os
from pathlib import Path
from typing import Literal, Optional
import numpy as np
import pyvista as pv
import open3d as o3d
from ribctl import RIBXZ_TEMP_FILES
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.alphalib import cif_to_point_cloud, fast_normal_estimation, quick_surface_points, validate_mesh_pyvista
from ribctl.lib.npet.kdtree_approach import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    apply_poisson_reconstruction,
    clip_tunnel_by_chain_proximity,
    create_point_cloud_mask,
    estimate_normals,
    filter_residues_parallel,
    ptcloud_convex_hull_points,
    ribosome_entities,
    transform_points_from_C0,
    transform_points_to_C0,
)
from ribctl.lib.npet.pipeline.pipeline_status_tracker import NPETProcessingTracker, ProcessingStage
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager
from ribctl.lib.npet.various_visualization import (
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
    visualize_filtered_residues,
    visualize_mesh,
    visualize_pointcloud,
    visualize_pointcloud_axis,
)
from ribctl.lib.schema.types_ribosome import ConstrictionSite, PTCInfo
from ribctl.ribosome_ops import RibosomeOps


def create_npet_mesh(rcsb_id: str, log_dir: Optional[Path] = None, force: bool = False) -> NPETProcessingTracker:
    """
    Creates NPET mesh for a given RCSB ID with comprehensive tracking and logging.
    
    Args:
        rcsb_id: The RCSB PDB identifier
        log_dir: Directory to store logs (default: logs)
        force: Whether to force regeneration of existing assets
        
    Returns:
        NPETProcessingTracker: Processing tracker with detailed status information
    """
    print(f"Creating NPET mesh for {rcsb_id}")
    rcsb_id = rcsb_id.upper()
    tracker = NPETProcessingTracker(rcsb_id, log_dir)
    
    try:
        # SETUP
        tracker.begin_stage(ProcessingStage.SETUP)
        assets = TunnelMeshAssetsManager(rcsb_id)
        cifpath = AssetType.MMCIF.get_path(rcsb_id)
        ashapepath = AssetType.ALPHA_SHAPE.get_path(rcsb_id)
        meshpath = AssetType.NPET_MESH.get_path(rcsb_id)

        ro = RibosomeOps(rcsb_id)
        profile = ro.profile
        
        # Create artifacts directory for this structure if it doesn't exist
        structure_dir = Path(ro.assets.paths.dir)
        artifacts_dir = structure_dir / "artifacts"
        artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        # Core parameters that affect multiple stages
        R = 35
        H = 120
        Vsize = 1
        ATOM_SIZE = 2

        normals_pcd_path = assets.tunnel_pcd_normal_estimated

        tracker.add_artifact(ProcessingStage.SETUP, Path(cifpath))
        tracker.end_stage(ProcessingStage.SETUP, True)
        
        # PTC IDENTIFICATION
        tracker.begin_stage(ProcessingStage.PTC_IDENTIFICATION)
        try:
            ptc_info = PTC_location(rcsb_id)
            ptc_pt = np.array(ptc_info.location)
            
            # Save PTC info as JSON
            ptc_json_path = AssetType.PTC.get_path(rcsb_id)
            with open(ptc_json_path, 'w') as f:
                f.write(ptc_info.model_dump_json())
            
            tracker.add_artifact(ProcessingStage.PTC_IDENTIFICATION, ptc_json_path)
            tracker.end_stage(ProcessingStage.PTC_IDENTIFICATION, True)
        except Exception as e:
            tracker.end_stage(ProcessingStage.PTC_IDENTIFICATION, False, error=e)
            tracker.complete_processing(False)
            return tracker
            
        # CONSTRICTION IDENTIFICATION
        tracker.begin_stage(ProcessingStage.CONSTRICTION_IDENTIFICATION)
        try:
            constriction_pt = get_constriction(rcsb_id)
            
            # Save constriction point as numpy array
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(rcsb_id)
            with open(constriction_path, 'w') as f:
                json.dump(ConstrictionSite(location=constriction_pt.tolist()).model_dump(), f)
            
            tracker.add_artifact(ProcessingStage.CONSTRICTION_IDENTIFICATION, constriction_path)
            tracker.end_stage(ProcessingStage.CONSTRICTION_IDENTIFICATION, True)
        except Exception as e:
            tracker.end_stage(ProcessingStage.CONSTRICTION_IDENTIFICATION, False, error=e)
            tracker.complete_processing(False)
            return tracker
        
        # ALPHA SHAPE GENERATION
        alpha_params = {
            "d3d_alpha": 75,    
            "d3d_tol": 4,      
            "d3d_offset": 3,    
            "kdtree_radius": 30,    
            "max_nn": 40,          
            "tangent_planes_k": 15,
            "PR_depth": 6,       
            "PR_ptweight": 4,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.ALPHA_SHAPE, alpha_params)
        
        if should_process or force:
            try:
                # Check if alpha shape already exists and should be used
                regenerate_alpha = force or not os.path.exists(ashapepath)
                
                if not regenerate_alpha:
                    print(f"Using existing alpha shape mesh: {ashapepath}")
                    tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, Path(ashapepath))
                else:
                    print(f"Generating alpha shape mesh for {rcsb_id}")
                    
                    # Generate the point cloud for the structure
                    ptcloudpath = artifacts_dir / f"{rcsb_id}_structure_ptcloud.npy"
                    
                    if not os.path.exists(ptcloudpath) or force:
                        print("Extracting point cloud from CIF file")
                        first_assembly_chains = ro.first_assembly_auth_asym_ids()
                        ptcloud = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
                        np.save(ptcloudpath, ptcloud)
                    else:
                        print(f"Loading existing point cloud from {ptcloudpath}")
                        ptcloud = np.load(ptcloudpath)
                    
                    tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, ptcloudpath)
                    
                    # Parameters for alpha shape generation
                    d3d_alpha = alpha_params["d3d_alpha"]
                    d3d_tol = alpha_params["d3d_tol"]
                    d3d_offset = alpha_params["d3d_offset"]
                    
                    kdtree_radius = alpha_params["kdtree_radius"]
                    max_nn = alpha_params["max_nn"]
                    tangent_planes_k = alpha_params["tangent_planes_k"]
                    
                    PR_depth = alpha_params["PR_depth"]
                    PR_ptweight = alpha_params["PR_ptweight"]
                    
                    print("Beginning Delaunay 3D reconstruction")
                    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)
                    
                    # Save surface points
                    surface_pts_path = artifacts_dir / f"{rcsb_id}_alpha_surface_points.npy"
                    np.save(surface_pts_path, surface_pts)
                    tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, surface_pts_path)
                    
                    # Visualize if possible
                    try:
                        surface_viz_path = artifacts_dir / f"{rcsb_id}_alpha_surface.png"
                        visualize_pointcloud(surface_pts, rcsb_id, output_path=str(surface_viz_path))
                        tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, surface_viz_path)
                    except Exception as viz_error:
                        print(f"Warning: Could not save alpha surface visualization: {str(viz_error)}")
                    
                    # Normal estimation
                    normal_estimated_pcd = fast_normal_estimation(
                        surface_pts, kdtree_radius, max_nn, tangent_planes_k
                    )
                    
                    # Save normal-estimated point cloud
                    alpha_normals_path = artifacts_dir / f"{rcsb_id}_alpha_normals.ply"
                    o3d.io.write_point_cloud(str(alpha_normals_path), normal_estimated_pcd)
                    tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, alpha_normals_path)
                    
                    # Apply Poisson reconstruction
                    apply_poisson_reconstruction(
                        str(alpha_normals_path),
                        ashapepath,
                        recon_depth=PR_depth,
                        recon_pt_weight=PR_ptweight,
                    )
                    
                    # Extract largest component for better quality
                    mesh = pv.read(ashapepath)
                    labeled = mesh.connectivity(largest=True)  # This keeps only the largest component
                    labeled.save(ashapepath)
                    
                    # Check watertightness
                    watertight = validate_mesh_pyvista(labeled)
                    
                    if not watertight:
                        print("Warning: Alpha shape mesh is not watertight, but continuing anyway")
                    
                    # Save mesh visualization
                    try:
                        mesh_viz_path = artifacts_dir / f"{rcsb_id}_alpha_mesh.png"
                        visualize_mesh(ashapepath, rcsb_id, output_path=str(mesh_viz_path))
                        tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, mesh_viz_path)
                    except Exception as viz_error:
                        print(f"Warning: Could not save alpha mesh visualization: {str(viz_error)}")
                    
                    # Also add ASCII version as artifact if it exists
                    ascii_path = Path(str(ashapepath).split(".")[0] + "_ascii.ply")
                    if ascii_path.exists():
                        tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, ascii_path)
                
                # Add the final alpha shape as artifact
                tracker.add_artifact(ProcessingStage.ALPHA_SHAPE, Path(ashapepath))
                tracker.end_stage(ProcessingStage.ALPHA_SHAPE, True)
                
                # Verify alpha shape exists for the remaining pipeline
                if not os.path.exists(ashapepath):
                    raise FileNotFoundError(f"Alpha shape file {ashapepath} not found after generation")
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.ALPHA_SHAPE, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # LANDMARK IDENTIFICATION
        # This stage is fast, so we don't need parameter tracking
        tracker.begin_stage(ProcessingStage.LANDMARK_IDENTIFICATION)
        try:
            # Load PTC and constriction points from previous stages
            ptc_path = AssetType.PTC.get_path(rcsb_id)
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(rcsb_id)
            
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
    
        except Exception as e:
            tracker.end_stage(ProcessingStage.LANDMARK_IDENTIFICATION, False, error=e)
            tracker.complete_processing(False)
            return tracker
            
        # ENTITY FILTERING
        # Define parameters for entity filtering
        entity_filter_params = {
            "radius": R,
            "height": H,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.ENTITY_FILTERING, entity_filter_params)
        
        if should_process or force:
            try:
                tunnel_debris = {
                    "3J7Z": ["a", "7"],
                    "7A5G": ["Y2"],
                    "5GAK": ["z"],
                    "5NWY": ["s"],
                }
                
                if profile.mitochondrial:
                    try:
                        chain = ro.get_poly_by_polyclass("mL45")
                        tunnel_debris[rcsb_id] = [chain.auth_asym_id]    
                    except Exception as chain_error:
                        print("Mitochondrial mL45 chain not found:", str(chain_error))
                
                residues = ribosome_entities(
                    rcsb_id,
                    cifpath,
                    "R",
                    tunnel_debris[rcsb_id] if rcsb_id in tunnel_debris else [],
                )
                
                filtered_residues = filter_residues_parallel(
                    residues, ptc_pt, constriction_pt, R, H
                )
                
                filtered_points = np.array(
                    [
                        atom.get_coord()
                        for residue in filtered_residues
                        for atom in residue.child_list
                    ]
                )
                
                # Save filtered points
                filtered_points_path = artifacts_dir / f"{rcsb_id}_filtered_points.npy"
                np.save(filtered_points_path, filtered_points)
                tracker.add_artifact(ProcessingStage.ENTITY_FILTERING, filtered_points_path)
                
                # Save visualization if possible
                try:
                    viz_path = artifacts_dir / f"{rcsb_id}_filtered_residues.png"
                    visualize_filtered_residues(
                        filtered_residues, residues, ptc_pt, constriction_pt, R, H,
                        output_path=str(viz_path)
                    )
                    tracker.add_artifact(ProcessingStage.ENTITY_FILTERING, viz_path)
                except Exception as viz_error:
                    print(f"Warning: Could not save visualization: {str(viz_error)}")
                
                tracker.end_stage(ProcessingStage.ENTITY_FILTERING, True)
            
            except Exception as e:
                tracker.end_stage(ProcessingStage.ENTITY_FILTERING, False, error=e)
                tracker.complete_processing(False)
                return tracker
            
        # POINT CLOUD PROCESSING
        ptcloud_params = {
            "radius": R,
            "height": H,
            "voxel_size": Vsize,
            "atom_size": ATOM_SIZE,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.POINT_CLOUD_PROCESSING, ptcloud_params)
        
        if should_process or force:
            try:
                transformed_points = transform_points_to_C0(
                    filtered_points, ptc_pt, constriction_pt
                )
                
                # Save transformed points
                transformed_points_path = artifacts_dir / f"{rcsb_id}_transformed_points.npy"
                np.save(transformed_points_path, transformed_points)
                tracker.add_artifact(ProcessingStage.POINT_CLOUD_PROCESSING, transformed_points_path)
                
                mask, (x, y, z) = create_point_cloud_mask(
                    transformed_points,
                    radius=R,
                    height=H,
                    voxel_size=Vsize,
                    radius_around_point=ATOM_SIZE,
                )
                
                points = np.where(~mask)
                empty_coordinates = np.column_stack((x[points[0]], y[points[1]], z[points[2]]))
                back_projected = transform_points_from_C0(
                    empty_coordinates, ptc_pt, constriction_pt
                )
                
                # Save back projected points
                back_projected_path = artifacts_dir / f"{rcsb_id}_back_projected.npy"
                np.save(back_projected_path, back_projected)
                tracker.add_artifact(ProcessingStage.POINT_CLOUD_PROCESSING, back_projected_path)
                
                ashape_watertight_mesh = pv.read(ashapepath)
                select = pv.PolyData(back_projected).select_enclosed_points(ashape_watertight_mesh)
                mask = select["SelectedPoints"]
                interior = back_projected[mask == 1]
                empty_in_world_coords = np.array(interior)
                
                # Save point cloud
                interior_points_path = artifacts_dir / f"{rcsb_id}_interior_points.npy"
                np.save(interior_points_path, empty_in_world_coords)
                tracker.add_artifact(ProcessingStage.POINT_CLOUD_PROCESSING, interior_points_path)
                
                # Save visualization if possible
                try:
                    pc_viz_path = artifacts_dir / f"{rcsb_id}_point_cloud.png"
                    visualize_pointcloud(empty_in_world_coords, rcsb_id, output_path=str(pc_viz_path))
                    tracker.add_artifact(ProcessingStage.POINT_CLOUD_PROCESSING, pc_viz_path)
                except Exception as viz_error:
                    print(f"Warning: Could not save point cloud visualization: {str(viz_error)}")
                
                tracker.end_stage(ProcessingStage.POINT_CLOUD_PROCESSING, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.POINT_CLOUD_PROCESSING, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # CLUSTERING
        clustering_params = {
            "epsilon": 5.5,
            "min_samples": 600,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.CLUSTERING, clustering_params)
        
        if should_process or force:
            try:
                _u_EPSILON_initial_pass = clustering_params["epsilon"]
                _u_MIN_SAMPLES_initial_pass = clustering_params["min_samples"]
                
                db, clusters_container = DBSCAN_capture(
                    empty_in_world_coords, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass
                )
                
                largest_cluster, largest_cluster_id = DBSCAN_pick_largest_cluster(
                    clusters_container
                )
                
                # Save largest cluster
                largest_cluster_path = artifacts_dir / f"{rcsb_id}_largest_cluster.npy"
                np.save(largest_cluster_path, largest_cluster)
                tracker.add_artifact(ProcessingStage.CLUSTERING, largest_cluster_path)
                
                # Save visualization if possible
                try:
                    cluster_viz_path = artifacts_dir / f"{rcsb_id}_clusters.png"
                    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                        clusters_container,
                        _u_EPSILON_initial_pass,
                        _u_MIN_SAMPLES_initial_pass,
                        ptc_pt,
                        constriction_pt,
                        largest_cluster,
                        R,
                        H,
                        output_path=str(cluster_viz_path)
                    )
                    tracker.add_artifact(ProcessingStage.CLUSTERING, cluster_viz_path)
                except Exception as viz_error:
                    print(f"Warning: Could not save cluster visualization: {str(viz_error)}")
                
                # Save clusters data
                try:
                    # Save each cluster as a separate file
                    for cluster_id, cluster_points in clusters_container.items():
                        if cluster_id != -1:  # Skip noise
                            cluster_path = artifacts_dir / f"{rcsb_id}_cluster_{cluster_id}.npy"
                            np.save(cluster_path, np.array(cluster_points))
                            if cluster_id == largest_cluster_id:
                                tracker.add_artifact(ProcessingStage.CLUSTERING, cluster_path)
                except Exception as cluster_error:
                    print(f"Warning: Could not save individual clusters: {str(cluster_error)}")
                
                tracker.end_stage(ProcessingStage.CLUSTERING, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.CLUSTERING, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # REFINEMENT
        refinement_params = {
            "epsilon": 3.5,
            "min_samples": 175,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.REFINEMENT, refinement_params)
        
        if should_process or force:
            try:
                _u_EPSILON_refinement = refinement_params["epsilon"]
                _u_MIN_SAMPLES_refinement = refinement_params["min_samples"]
                
                db_2, refined_clusters_container = DBSCAN_capture(
                    largest_cluster, _u_EPSILON_refinement, _u_MIN_SAMPLES_refinement
                )
                
                refined_cluster, refined_cluster_id = DBSCAN_pick_largest_cluster(
                    refined_clusters_container
                )
                
                # Save refined cluster
                refined_cluster_path = artifacts_dir / f"{rcsb_id}_refined_cluster.npy"
                np.save(refined_cluster_path, refined_cluster)
                tracker.add_artifact(ProcessingStage.REFINEMENT, refined_cluster_path)
                
                # Save visualization if possible
                try:
                    refined_viz_path = artifacts_dir / f"{rcsb_id}_refined_clusters.png"
                    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                        refined_clusters_container,
                        _u_EPSILON_refinement,
                        _u_MIN_SAMPLES_refinement,
                        ptc_pt,
                        constriction_pt,
                        refined_cluster,
                        R,
                        H,
                        output_path=str(refined_viz_path)
                    )
                    tracker.add_artifact(ProcessingStage.REFINEMENT, refined_viz_path)
                except Exception as viz_error:
                    print(f"Warning: Could not save refined cluster visualization: {str(viz_error)}")
                
                tracker.end_stage(ProcessingStage.REFINEMENT, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.REFINEMENT, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # SURFACE EXTRACTION
        surface_params = {
            "alpha": 2,
            "tolerance": 1,
            "offset": 2,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.SURFACE_EXTRACTION, surface_params)
        
        if should_process or force:
            try:
                d3d_alpha = surface_params["alpha"]
                d3d_tol = surface_params["tolerance"]
                d3d_offset = surface_params["offset"]
                
                surface_pts = ptcloud_convex_hull_points(
                    refined_cluster, d3d_alpha, d3d_tol, d3d_offset
                )
                
                # Save surface points
                surface_pts_path = artifacts_dir / f"{rcsb_id}_surface_points.npy"
                np.save(surface_pts_path, surface_pts)
                tracker.add_artifact(ProcessingStage.SURFACE_EXTRACTION, surface_pts_path)
                
                # Save visualization if possible
                try:
                    surface_viz_path = artifacts_dir / f"{rcsb_id}_surface_points.png"
                    visualize_pointcloud(surface_pts, rcsb_id, output_path=str(surface_viz_path))
                    tracker.add_artifact(ProcessingStage.SURFACE_EXTRACTION, surface_viz_path)
                except Exception as viz_error:
                    print(f"Warning: Could not save surface point visualization: {str(viz_error)}")
                
                tracker.end_stage(ProcessingStage.SURFACE_EXTRACTION, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.SURFACE_EXTRACTION, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # NORMAL ESTIMATION
        normal_params = {
            "kdtree_radius": 10,
            "kdtree_max_nn": 15,
            "correction_tangent_planes_n": 10,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.NORMAL_ESTIMATION, normal_params)
        
        if should_process or force:
            try:
                normal_estimated_pcd = estimate_normals(
                    surface_pts,
                    kdtree_radius=normal_params["kdtree_radius"],
                    kdtree_max_nn=normal_params["kdtree_max_nn"],
                    correction_tangent_planes_n=normal_params["correction_tangent_planes_n"],
                )
                
                # Instead of using a temp file, save in the structure's directory
                normals_pcd_path = artifacts_dir / f"{rcsb_id}_normal_estimated_pcd.ply"
                o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
                tracker.add_artifact(ProcessingStage.NORMAL_ESTIMATION, normals_pcd_path)
                
                tracker.end_stage(ProcessingStage.NORMAL_ESTIMATION, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.NORMAL_ESTIMATION, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # MESH RECONSTRUCTION
        mesh_params = {
            "depth": 6,
            "ptweight": 3,
        }
        
        should_process = tracker.begin_stage(ProcessingStage.MESH_RECONSTRUCTION, mesh_params)
        
        if should_process or force:
            try:
                PR_depth = mesh_params["depth"]
                PR_ptweight = mesh_params["ptweight"]
                
                apply_poisson_reconstruction(
                    str(normals_pcd_path),  # Use the normals path from artifacts_dir
                    meshpath,
                    recon_depth=PR_depth,
                    recon_pt_weight=PR_ptweight,
                )
                
                # Add mesh as artifact
                tracker.add_artifact(ProcessingStage.MESH_RECONSTRUCTION, Path(meshpath))
                
                # Also save ASCII version as artifact if it exists
                ascii_path = Path(str(meshpath).split(".")[0] + "_ascii.ply")
                if ascii_path.exists():
                    tracker.add_artifact(ProcessingStage.MESH_RECONSTRUCTION, ascii_path)
                
                tracker.end_stage(ProcessingStage.MESH_RECONSTRUCTION, True)
                
            except Exception as e:
                tracker.end_stage(ProcessingStage.MESH_RECONSTRUCTION, False, error=e)
                tracker.complete_processing(False)
                return tracker
        
        # VALIDATION
        tracker.begin_stage(ProcessingStage.VALIDATION)
        try:
            # Save mesh visualization if possible
            try:
                mesh_viz_path = artifacts_dir / f"{rcsb_id}_mesh.png"
                visualize_mesh(meshpath, rcsb_id, output_path=str(mesh_viz_path))
                tracker.add_artifact(ProcessingStage.VALIDATION, mesh_viz_path)
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
            
            tracker.end_stage(ProcessingStage.VALIDATION, True)
            tracker.complete_processing(True, watertight=True)
            
        except Exception as e:
            tracker.end_stage(ProcessingStage.VALIDATION, False, error=e)
            tracker.complete_processing(False, watertight=False)
            return tracker
        
        return tracker
            
    except Exception as e:
        # Catch any unhandled exceptions
        if tracker.current_stage:
            tracker.end_stage(tracker.current_stage, False, error=e)
        tracker.complete_processing(False)
        return tracker
