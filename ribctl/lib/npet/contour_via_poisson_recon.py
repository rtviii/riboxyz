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
    """
    Convert a CIF file to a point cloud, optionally filtering for specific chains.

    Args:
        cif_path (str): Path to the CIF file
        chains (list[str] | None): Optional list of chain IDs to include. If None, includes all chains.

    """
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




def alpha_contour_via_poisson_recon(rcsb_id:str, verbose:bool=False):
    rops                  = RibosomeOps(rcsb_id)
    cifpath               = rops.assets.paths.cif
    first_assembly_chains = rops.first_assembly_auth_asym_ids()
    ptcloud               = cif_to_point_cloud(cifpath, first_assembly_chains, do_atoms=True)

    if verbose:
        visualize_pointcloud(ptcloud)

    output_normals_pcd = os.path.join(RIBXZ_TEMP_FILES, "{}_normal_estimated_pcd.ply")
    output_mesh        = AssetType.ALPHA_SHAPE.get_path(rcsb_id)

    d3d_alpha  = 15  # Increased from 2
    d3d_tol    = 1
    d3d_offset = 2

    kdtree_radius   = 10    # Reduced from 10
    max_nn          = 20         # Reduced from 100
    tanget_planes_k = 10

    PR_depth    = 6        # Reduced from 6
    PR_ptweight = 4

    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)
    if verbose:
        visualize_pointcloud(surface_pts, rcsb_id)

    normal_estimated_pcd = fast_normal_estimation(surface_pts, kdtree_radius, max_nn, tanget_planes_k)

    if verbose:
        o3d.visualization.draw_geometries([normal_estimated_pcd], point_show_normal=True)

    o3d.io.write_point_cloud(output_normals_pcd, normal_estimated_pcd)

    apply_poisson_reconstruction(
        output_normals_pcd,
        Path(output_mesh),
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    validate_mesh_pyvista(pv.read(output_mesh))
    if verbose:
        visualize_mesh(output_mesh)
