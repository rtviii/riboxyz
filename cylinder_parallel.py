import os
import numpy as np
from typing import Tuple
import multiprocessing as mp
from cylinder import transform_points_to_C0
from mesh_generation.mes_visualization import visualize_pointcloud
from numba import jit
import pyvista as pv

from ribctl.lib.landmarks.constriction import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import filter_residues_parallel, ribosome_entities

def chunk_points(points: np.ndarray, n_chunks: int) -> list:
    """Split points into chunks for parallel processing"""
    chunk_size = len(points) // n_chunks
    return [points[i:i + chunk_size] for i in range(0, len(points), chunk_size)]

@jit(nopython=True)
def process_points_chunk(points: np.ndarray, grid_coords: np.ndarray, radius: float) -> np.ndarray:
    """Process a chunk of points with Numba acceleration"""
    mask = np.zeros(grid_coords.shape[1], dtype=np.bool_)
    for point in points:
        distances = np.sqrt(np.sum((grid_coords.T - point)**2, axis=1))
        mask |= (distances <= radius)
    return mask

def parallel_point_cloud_mask(points: np.ndarray, X: np.ndarray, Y: np.ndarray, Z: np.ndarray, 
                            radius_around_point: float) -> np.ndarray:
    """Generate point cloud mask using parallel processing"""
    # Prepare grid coordinates once
    grid_coords = np.stack([X, Y, Z])
    original_shape = X.shape
    grid_coords = grid_coords.reshape(3, -1)
    
    # Determine number of chunks based on CPU cores
    n_cores = mp.cpu_count() - 4 
    chunks = chunk_points(points, n_cores)
    
    # Process chunks in parallel
    with mp.Pool(n_cores) as pool:
        results = pool.starmap(process_points_chunk, 
                             [(chunk, grid_coords, radius_around_point) for chunk in chunks])
    
    # Combine results
    final_mask = np.any(results, axis=0)
    return final_mask.reshape(original_shape)

def main():
    # Your existing setup code...
    RCSB_ID = '3J7Z'
    radius     = 40
    height     = 80
    voxel_size = 1
    ATOM_RADIUS = 2

    base_point = np.array(PTC_location(RCSB_ID).location)
    axis_point = np.array( get_constriction(RCSB_ID) )

    if os.path.exists('points.npy'):
        points = np.load('points.npy')
        print("Loaded")
    else:
        residues= filter_residues_parallel( ribosome_entities(RCSB_ID, 'R'), base_point, axis_point, radius, height, )
        points = np.array([atom.get_coord() for residue in residues for atom in residue.child_list])
        np.save('points.npy', points)
        print("Saved")

    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

    # Transform points (vectorized)
    transformed = transform_points_to_C0(points, base_point, axis_point)
    X_I, Y_I, Z_I = transformed.T
    points = np.column_stack((X_I, Y_I, Z_I))

    # Generate cylinder mask (vectorized)
    cylinder_mask = (np.sqrt(X**2 + Y**2) <= radius)
    hollow_cylinder = ~cylinder_mask

    # Generate point cloud mask in parallel
    point_cloud_mask = parallel_point_cloud_mask(points, X, Y, Z, ATOM_RADIUS)

    # Combine masks
    final_mask = hollow_cylinder | point_cloud_mask

    # Visualize results
    occupied = np.where(~final_mask)
    visualization_points = np.column_stack((
        x[occupied[0]], 
        y[occupied[1]], 
        z[occupied[2]]
    ))
    occupied_points = pv.PolyData(visualization_points)
    visualize_pointcloud(occupied_points)


if __name__ == '__main__':
    main()