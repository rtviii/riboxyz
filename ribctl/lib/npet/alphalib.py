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

def cif_to_point_cloud(cif_path: str, chains: list[str] | None = None,  do_atoms:bool=False, exclude_chains: list[str] = []) -> np.ndarray:
    print("Excluding chains:", exclude_chains)
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    coordinates = []

    first_model = structure[0]
    if do_atoms:
        for chain in first_model:
            if exclude_chains is not None and chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.get_coord())    
    else:
        for chain in first_model:
            if exclude_chains is not None and chain.id in exclude_chains:
                continue
            if chains is not None and chain.id not in chains:
                continue
            for residue in chain:
                coordinates.append(residue.center_of_mass())

    if not coordinates:
        raise ValueError(f"No coordinates found in {cif_path}")

    return np.array(coordinates)

def validate_mesh_pyvista(mesh_or_path, stage="unknown"):
    """
    Validates if a mesh is watertight.
    
    Parameters:
        mesh_or_path: Either a PyVista mesh object or a path to a mesh file
        stage: Optional identifier for logging
        
    Returns:
        bool: True if the mesh is watertight, False otherwise
    """
    import pyvista as pv
    import os
    from pathlib import Path

    # If mesh is a path, load it
    if isinstance(mesh_or_path, (str, Path)):
        if not os.path.exists(mesh_or_path):
            print(f"WARNING: Mesh file does not exist: {mesh_or_path}")
            return False
        
        try:
            mesh = pv.read(str(mesh_or_path))
        except Exception as e:
            print(f"WARNING: Failed to load mesh at {mesh_or_path}: {e}")
            return False
    else:
        mesh = mesh_or_path

    if mesh is None:
        print(f"WARNING: Null mesh at stage {stage}")
        return False

    print(f"\nMesh properties at stage: {stage}")

    try:
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
    except Exception as e:
        print(f"WARNING: Error checking watertightness: {e}")
        return False

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


    output_normals_pcd = os.path.join(RIBXZ_TEMP_FILES, "{}_normal_estimated_pcd.ply".format(rcsb_id))
    output_mesh        = AssetType.ALPHA_SHAPE.get_path(rcsb_id)

    d3d_alpha  = 75    # Increase from 35 - be more aggressive
    d3d_tol    = 4     # Increase tolerance to smooth out small details
    d3d_offset = 3     # Slightly larger offset

    kdtree_radius   = 30    # Larger radius to catch more global structure
    max_nn          = 40    # More neighbors for more robust estimation
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