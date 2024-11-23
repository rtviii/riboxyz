import numpy as np
from scipy.spatial import cKDTree
import pyvista as pv

from cylinder import get_transformation_to_C0, transform_points_to_C0
from mesh_generation.mes_visualization import visualize_pointcloud
from ribctl.lib.landmarks.constriction import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import filter_residues_parallel, ribosome_entities

def generate_voxel_centers(radius: float, height: float, voxel_size: float) -> tuple:
    """Generate centers of all voxels in the grid"""
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)
    
    # Generate all voxel center coordinates
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    
    return voxel_centers, (X.shape, x, y, z)

def create_point_cloud_mask(points: np.ndarray, 
                          radius: float, 
                          height: float,
                          voxel_size: float = 1.0,
                          radius_around_point: float = 2.0):
    """
    Create point cloud mask using KDTree for efficient spatial queries
    """
    # Generate voxel centers
    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(radius, height, voxel_size)
    
    # Create KDTree from the transformed points
    tree = cKDTree(points)
    
    # Find all voxels that have points within radius_around_point
    # This is much more efficient than checking each point against each voxel
    indices = tree.query_ball_point(voxel_centers, radius_around_point)
    
    # Create mask from the indices
    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True
    
    # Reshape mask back to grid shape
    point_cloud_mask = point_cloud_mask.reshape(grid_shape)
    
    # Create cylinder mask
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    cylinder_mask = (np.sqrt(X**2 + Y**2) <= radius)
    hollow_cylinder = ~cylinder_mask
    
    # Combine masks
    final_mask = hollow_cylinder | point_cloud_mask
    
    return final_mask, (x, y, z)

def transform_points_to_C0(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_translated  = points + translation
    points_transformed = points_translated @ rotation.T
    
    return points_transformed


def transform_points_from_C0(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_unrotated = points @ rotation
    points_untranslated = points_unrotated - translation
    
    return points_untranslated
def main():
    # Load your points and transform them as before
    RCSB_ID    = '4UG0'
    R          = 15
    H          = 120
    Vsize      = 1
    ATOM_SIZE  = 2
    base_point = np.array(PTC_location(RCSB_ID).location)
    axis_point = np.array(get_constriction(RCSB_ID) )

    residues           = filter_residues_parallel( ribosome_entities(RCSB_ID, 'R'), base_point, axis_point, R, H)
    points             = np.array([atom.get_coord() for residue in residues for atom in residue.child_list])
    transformed_points = transform_points_to_C0(points, base_point, axis_point)

    mask, (x, y, z) = create_point_cloud_mask(
        transformed_points,
        radius              = R,
        height              = H,
        voxel_size          = Vsize,
        radius_around_point = ATOM_SIZE
    )
    
    points = np.where(~mask)
    empty_coordinates = np.column_stack((
        x[points[0]], 
        y[points[1]], 
        z[points[2]]
    ))
    world_coords    = transform_points_from_C0(empty_coordinates ,base_point,axis_point)
    occupied_points = pv.PolyData(empty_coordinates)
    world_coords    = pv.PolyData(world_coords)

    visualize_pointcloud(occupied_points, world_coords)

if __name__ == '__main__':
    main()