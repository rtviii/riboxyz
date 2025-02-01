from pathlib import Path
import sys

sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.lib.npet.kdtree_approach import apply_poisson_reconstruction

from Bio.PDB.MMCIFParser import MMCIFParser
from alpha_lib import cif_to_point_cloud, validate_mesh_pyvista
from ribctl.lib.npet.various_visualization import visualize_mesh, visualize_pointcloud
from ribctl.ribosome_ops import RibosomeOps

import numpy as np
import pyvista as pv
import open3d as o3d


def validate_mesh_pyvista(mesh, stage="unknown"):
    """Validate and print mesh properties, focusing on watertightness."""
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
    return mesh


def quick_surface_points(
    pointcloud: np.ndarray, alpha: float, tolerance: float, offset: float
) -> np.ndarray:
    cloud = pv.PolyData(pointcloud)
    # Using larger tolerance and smaller offset for faster computation
    grid = cloud.delaunay_3d( alpha=alpha, tol=tolerance, offset=offset, progress_bar=True )
    surface = grid.extract_surface().cast_to_pointset()
    return surface.points


def fast_normal_estimation(
    surface_pts: np.ndarray,
    kdtree_radius,
    max_nn,
    tangent_planes_k, 
) -> o3d.geometry.PointCloud:
    """
    Estimate normals for surface points with optimized parameters for speed.

    Args:
        surface_pts: Input surface points
        kdtree_radius: Search radius for neighbors (smaller = faster)
        max_nn: Maximum number of neighbors to consider (smaller = faster)
        tangent_planes_k: Number of neighbors for tangent plane estimation
    Returns:
        o3d.geometry.PointCloud: Point cloud with estimated normals
    """
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(surface_pts)

    # Use hybrid search with reduced parameters
    search_param = o3d.geometry.KDTreeSearchParamHybrid(
        radius=kdtree_radius, max_nn=max_nn
    )

    # Estimate normals
    pcd.estimate_normals(search_param=search_param)

    # Orient normals with fewer neighbors
    pcd.orient_normals_consistent_tangent_plane(k=tangent_planes_k)

    return pcd


def create_mesh_from_pointcloud(
    pointcloud: np.ndarray,
    poisson_depth: int = 6,  # Reduced from typical 8-9
    poisson_weight: float = 2.0,
    visualize: bool = False,
) -> tuple[np.ndarray, o3d.geometry.PointCloud]:
    """
    Complete pipeline to create a mesh from a point cloud, optimized for speed.

    Args:
        pointcloud: Input point cloud
        poisson_depth: Depth for Poisson reconstruction (lower = faster)
        poisson_weight: Point weight for reconstruction
        visualize: Whether to visualize the normal estimation
    Returns:
        tuple: (surface points, point cloud with normals)
    """
    # Get surface points with relaxed parameters
    surface_pts = quick_surface_points(
        pointcloud,
        alpha=2.0,  # Increased alpha for faster computation
        tolerance=0.01,  # Increased tolerance for faster computation
    )

    # Estimate normals with reduced parameters
    pcd_with_normals = fast_normal_estimation(
        surface_pts,
        kdtree_radius=5.0,  # Reduced radius
        max_nn=10,  # Reduced neighbor count
        tangent_planes_k=5,  # Reduced tangent plane neighbors
    )

    if visualize:
        o3d.visualization.draw_geometries([pcd_with_normals], point_show_normal=True)

    return surface_pts, pcd_with_normals


if __name__ == "__main__":

    rops                  = RibosomeOps("4TUA")
    cifpath               = rops.assets.paths.cif
    first_assembly_chains = rops.first_assembly_auth_asym_ids()
    ptcloud               = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)
    # np.save("4TUA_pointcloud.npy", ptcloud)
    # ptcloud = np.load("4TUA_pointcloud.npy")
    visualize_pointcloud(ptcloud)

    output_normals_pcd = "output_normals_pcd.ply"
    output_mesh = "output_mesh.ply"

    d3d_alpha  = 15  # Increased from 2
    d3d_tol    = 1
    d3d_offset = 2

    # Normal estimation
    kdtree_radius   = 10    # Reduced from 10
    max_nn          = 20         # Reduced from 100
    tanget_planes_k = 10

    # Poisson reconstruction
    PR_depth = 6        # Reduced from 6
    PR_ptweight = 4   

    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)
    visualize_pointcloud(surface_pts, "4TUA")

    normal_estimated_pcd = fast_normal_estimation(surface_pts, kdtree_radius, max_nn, tanget_planes_k)
    o3d.visualization.draw_geometries([normal_estimated_pcd], point_show_normal=True)

    o3d.io.write_point_cloud(output_normals_pcd, normal_estimated_pcd)


    apply_poisson_reconstruction(
        output_normals_pcd,
        Path(output_mesh),
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )
    validate_mesh_pyvista(pv.read(output_mesh))
    visualize_mesh(output_mesh)
