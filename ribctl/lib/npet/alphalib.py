import os
from pathlib import Path
import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.asset_manager.asset_types import AssetType
from ribctl import RIBXZ_TEMP_FILES
from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction

from Bio.PDB.MMCIFParser import MMCIFParser
from ribctl.lib.npet.various_visualization import visualize_mesh, visualize_pointcloud
from ribctl.ribosome_ops import RibosomeOps

import numpy as np
import pyvista as pv
import open3d as o3d

def cif_to_point_cloud(cif_path: str, chains: list[str] | None = None,  do_atoms:bool=False):
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    coordinates = []

    first_model = structure[0]
    if do_atoms:
        for chain in first_model:
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())    
    else:
        for chain in first_model:
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                coordinates.append(residue.center_of_mass())

    if not coordinates:
        raise ValueError(f"No coordinates found in {cif_path}")

    return np.array(coordinates)

def validate_mesh_pyvista(mesh, stage="unknown"):
    if mesh is None:
        print(f"WARNING: Null mesh at stage {stage}")
        return None

    print(f"\nMesh properties at stage: {stage}")

    # Check watertightness by looking for boundary edges
    # A mesh is watertight if it has no boundary edges
    edges = mesh.extract_feature_edges(
        boundary_edges=True,
        feature_edges=False,
        manifold_edges=False,
        non_manifold_edges=False,
    )
    is_watertight = edges.n_cells == 0

    print(f"- Is watertight: {is_watertight}")
    print("DONE\n\n ")
    return is_watertight

def quick_surface_points(
    pointcloud: np.ndarray, alpha: float, tolerance: float, offset: float
) -> np.ndarray:
    cloud = pv.PolyData(pointcloud)
    # Using larger tolerance and smaller offset for faster computation
    grid    = cloud.delaunay_3d(alpha=alpha, tol=tolerance, offset=offset, progress_bar=True )
    surface = grid.extract_surface().cast_to_pointset()
    return surface.points

def fast_normal_estimation(
    surface_pts: np.ndarray,
    kdtree_radius,
    max_nn,
    tangent_planes_k, 
) -> o3d.geometry.PointCloud:
    """
    Returns:
        o3d.geometry.PointCloud: Point cloud with estimated normals
    """
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(surface_pts)

    # Use hybrid search with reduced parameters
    search_param = o3d.geometry.KDTreeSearchParamHybrid( radius=kdtree_radius, max_nn=max_nn )

    pcd.estimate_normals(search_param=search_param)
    pcd.orient_normals_consistent_tangent_plane(k=tangent_planes_k)

    return pcd

def alpha_contour_via_poisson_recon(rcsb_id:str, verbose:bool=False):
    ptcloudpath = os.path.join(RIBXZ_TEMP_FILES, '{}_ptcloud.npx'.format(rcsb_id))
    rops                  = RibosomeOps(rcsb_id)
    cifpath               = rops.assets.paths.cif


    print("Cifpath:", cifpath)
    if not os.path.exists(ptcloudpath):
        print("Extracting point cloud from CIF file")
        first_assembly_chains = rops.first_assembly_auth_asym_ids()
        ptcloud               = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
        np.save(ptcloudpath, ptcloud)
    else:
        print("Loaded.")
        ptcloud = np.load(ptcloudpath)

    print("Before visualize ptclodu")
    if verbose:
        visualize_pointcloud(ptcloud, rcsb_id)
    print("after visualize ptclodu")

    output_normals_pcd = os.path.join(RIBXZ_TEMP_FILES, "{}_normal_estimated_pcd.ply".format(rcsb_id))
    output_mesh        = AssetType.ALPHA_SHAPE.get_path(rcsb_id)

    d3d_alpha  = 45    # Increase from 35 - be more aggressive
    d3d_tol    = 2     # Increase tolerance to smooth out small details
    d3d_offset = 3     # Slightly larger offset

    kdtree_radius   = 30    # Larger radius to catch more global structure
    max_nn          = 30    # More neighbors for more robust estimation
    tanget_planes_k = 15    # Actually decrease this to avoid over-smoothing

    PR_depth    = 6        # Reduced from 6
    PR_ptweight = 4

    print("Beginning Delaunay 3d reconstruction")
    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)

    if verbose:
        visualize_pointcloud(surface_pts, rcsb_id)

    normal_estimated_pcd = fast_normal_estimation(surface_pts, kdtree_radius, max_nn, tanget_planes_k)

    if verbose:
        o3d.visualization.draw_geometries([normal_estimated_pcd], point_show_normal=True)

    o3d.io.write_point_cloud(output_normals_pcd, normal_estimated_pcd)

    apply_poisson_reconstruction(
        output_normals_pcd,
        output_mesh,
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    mesh = pv.read(output_mesh)
    labeled = mesh.connectivity(largest=True)  # This keeps only the largest component
    
    # Save the filtered mesh
    labeled.save(output_mesh)
    
    watertight = validate_mesh_pyvista(labeled)
    if verbose:
        visualize_mesh(output_mesh, rcsb_id)
    if not watertight:
        print("XXXX Watertightness check failed, removing", output_mesh , " XXXX")
        os.remove(output_mesh)



# import torch
# import kaolin
# from kaolin.ops.mesh import check_sign
# import numpy as np
# import pyvista as pv
# import kaolin.ops.conversions as cvt

# def compute_point_distances(query_points: torch.Tensor, reference_points: torch.Tensor) -> torch.Tensor:
#     """
#     Compute distances from query points to nearest reference points using GPU
#     """
#     # Compute pairwise distances
#     query_points = query_points.unsqueeze(1)      # Shape: (N, 1, 3)
#     reference_points = reference_points.unsqueeze(0)  # Shape: (1, M, 3)
    
#     # Efficient distance computation
#     distances = torch.sqrt(((query_points - reference_points) ** 2).sum(-1))
#     min_distances, _ = torch.min(distances, dim=1)
    
#     return min_distances

# def generate_watertight_mesh(
#     points: np.ndarray,
#     resolution: int = 64,
#     iso_value: float = 0.5,
#     hole_size: float = 100  # Maximum hole size to fill
# ) -> pv.PolyData:
#     """
#     Generate a watertight mesh using Kaolin's GPU implementation
    
#     Args:
#         points: Point cloud coordinates
#         resolution: Voxel grid resolution
#         iso_value: Surface threshold value
#         hole_size: Maximum hole size to fill
#     """
#     device = torch.device('cuda')
#     points_torch = torch.from_numpy(points).float().to(device)
    
#     # First convert points to voxel grid using Kaolin
#     voxels = cvt.pointclouds_to_voxelgrids(
#         points_torch.unsqueeze(0),  # Add batch dimension
#         resolution=resolution
#     )
    
#     # Convert voxels to mesh using Kaolin's marching cubes
#     vertices_list, faces_list = cvt.voxelgrids_to_trianglemeshes(
#         voxels,
#         iso_value=iso_value
#     )
    
#     # Get first (and only) mesh from batch
#     vertices = vertices_list[0].cpu().numpy()
#     faces = faces_list[0].cpu().numpy()
    
#     # Create PyVista mesh
#     mesh = pv.make_tri_mesh(vertices, faces)
#     # Clean and ensure watertightness
#     mesh = mesh.clean(tolerance=1e-6)
#     mesh = mesh.fill_holes(hole_size)  # Specify maximum hole size
#     mesh = mesh.triangulate()
#     mesh = mesh.clean()
    
#     return mesh


# def gpu_accelerated_watertight_mesh(rcsb_id: str, verbose: bool = False):
#     ptcloudpath = os.path.join(RIBXZ_TEMP_FILES, f'{rcsb_id}_ptcloud.npx')
#     rops = RibosomeOps(rcsb_id)
#     cifpath = rops.assets.paths.cif
    
#     if not os.path.exists(ptcloudpath):
#         first_assembly_chains = rops.first_assembly_auth_asym_ids()
#         ptcloud = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
#         np.save(ptcloudpath, ptcloud)
#         print("Generated and saved point cloud:", ptcloudpath)
#     else:
#         ptcloud = np.load(ptcloudpath)
    
#     output_mesh = AssetType.ALPHA_SHAPE.get_path(rcsb_id)
    
#     try:
#         mesh = generate_watertight_mesh(
#             ptcloud,
#             resolution = 32,
#             iso_value  = 0.75,
#             hole_size  = 5
#         )
        
#         # Save mesh
#         mesh.save(output_mesh)
        
#         # Validate
#         watertight = validate_mesh_pyvista(mesh)
#         if verbose:
#             visualize_mesh(output_mesh, rcsb_id)
        
#         if not watertight:
#             print("XXXX Watertightness check failed, removing", output_mesh, " XXXX")
#             os.remove(output_mesh)
#             return None
        
#         visualize_mesh(output_mesh, rcsb_id)
#         return mesh
        
#     except Exception as e:
#         print(f"Failed to generate mesh: {str(e)}")
#         if os.path.exists(output_mesh):
#             os.remove(output_mesh)
#         return None
    