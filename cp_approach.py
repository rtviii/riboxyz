import numpy as np
from scipy.spatial import cKDTree
import pyvista as pv
import cupy as cp
from numba import jit  # Only using Numba for CPU optimizations
from concurrent.futures import ThreadPoolExecutor
import math

from cylinder import get_transformation_to_C0
from ribctl.lib.landmarks.constriction import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import filter_residues_parallel, ribosome_entities

class OptimizedVoxelGrid:
    def __init__(self, use_gpu=True):
        self.use_gpu = use_gpu and cp.cuda.is_available()
        
    @staticmethod
    @jit(nopython=True, parallel=True)
    def _generate_voxel_centers_numba(radius: float, height: float, voxel_size: float):
        """Numba-accelerated voxel center generation"""
        nx = ny = int(2 * radius / voxel_size) + 1
        nz = int(height / voxel_size) + 1
        
        total_points = nx * ny * nz
        centers = np.empty((total_points, 3))
        
        for i in range(nx):
            x = -radius + i * voxel_size
            for j in range(ny):
                y = -radius + j * voxel_size
                for k in range(nz):
                    z = k * voxel_size
                    idx = i * ny * nz + j * nz + k
                    centers[idx] = np.array([x, y, z])
                    
        return centers, (nx, ny, nz)

    def _chunk_points(self, points, chunk_size=1000000):
        """Split points into manageable chunks for processing"""
        return np.array_split(points, max(1, len(points) // chunk_size))

    def create_point_cloud_mask(self, points: np.ndarray, 
                              radius: float, 
                              height: float,
                              voxel_size: float = 1.0,
                              radius_around_point: float = 2.0,
                              max_chunk_size: int = 1000000):
        """
        Create point cloud mask using parallel processing and GPU acceleration
        """
        # Generate voxel centers using Numba
        voxel_centers, grid_shape = self._generate_voxel_centers_numba(radius, height, voxel_size)
        
        # Create KDTree from points
        tree = cKDTree(points)
        
        # Process voxel centers in parallel using ThreadPoolExecutor
        def process_chunk(chunk_centers):
            indices = tree.query_ball_point(chunk_centers, radius_around_point)
            mask = np.zeros(len(chunk_centers), dtype=bool)
            mask[[i for i, idx in enumerate(indices) if idx]] = True
            return mask
        
        # Split voxel centers into chunks for parallel processing
        center_chunks = self._chunk_points(voxel_centers, max_chunk_size)
        
        # Process chunks in parallel
        with ThreadPoolExecutor() as executor:
            chunk_masks = list(executor.map(process_chunk, center_chunks))
        
        # Combine chunk results
        point_cloud_mask = np.concatenate(chunk_masks)
        
        # Reshape mask back to grid shape
        point_cloud_mask = point_cloud_mask.reshape(grid_shape)
        
        # Create cylinder mask (GPU-accelerated if available)
        if self.use_gpu:
            x = cp.linspace(-radius, radius, grid_shape[0])
            y = cp.linspace(-radius, radius, grid_shape[1])
            z = cp.linspace(0, height, grid_shape[2])
            X, Y, Z = cp.meshgrid(x, y, z, indexing='ij')
            cylinder_mask = (cp.sqrt(X**2 + Y**2) <= radius)
            hollow_cylinder = ~cylinder_mask
            hollow_cylinder = cp.asnumpy(hollow_cylinder)
        else:
            x = np.linspace(-radius, radius, grid_shape[0])
            y = np.linspace(-radius, radius, grid_shape[1])
            z = np.linspace(0, height, grid_shape[2])
            X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
            cylinder_mask = (np.sqrt(X**2 + Y**2) <= radius)
            hollow_cylinder = ~cylinder_mask
        
        # Combine masks
        final_mask = hollow_cylinder | point_cloud_mask
        
        return final_mask, (x, y, z)

def transform_points_optimized(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray, 
                             use_gpu: bool = True) -> np.ndarray:
    """Optimized point transformation using CuPy for GPU acceleration"""
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    
    if use_gpu and cp.cuda.is_available():
        # Move data to GPU
        points_gpu = cp.asarray(points)
        translation_gpu = cp.asarray(translation)
        rotation_gpu = cp.asarray(rotation)
        
        # Perform transformation on GPU
        points_translated = points_gpu + translation_gpu
        points_transformed = cp.matmul(points_translated, rotation_gpu.T)
        
        # Move result back to CPU
        return cp.asnumpy(points_transformed)
    else:
        # CPU implementation using NumPy
        points_translated = points + translation
        points_transformed = points_translated @ rotation.T
        return points_transformed

def transform_points_from_C0_optimized(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray,
                                     use_gpu: bool = True) -> np.ndarray:
    """Optimized inverse point transformation using CuPy for GPU acceleration"""
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    
    if use_gpu and cp.cuda.is_available():
        # Move data to GPU
        points_gpu = cp.asarray(points)
        translation_gpu = cp.asarray(translation)
        rotation_gpu = cp.asarray(rotation)
        
        # Perform inverse transformation on GPU
        points_unrotated = cp.matmul(points_gpu, rotation_gpu)
        points_untranslated = points_unrotated - translation_gpu
        
        # Move result back to CPU
        return cp.asnumpy(points_untranslated)
    else:
        # CPU implementation using NumPy
        points_unrotated = points @ rotation
        points_untranslated = points_unrotated - translation
        return points_untranslated

def main():
    # Example usage with your existing parameters
    RCSB_ID = '4UG0'
    R = 40
    H = 120
    Vsize = 0.5  # Twice as fine as original
    ATOM_SIZE = 2
    base_point = np.array(PTC_location(RCSB_ID).location)
    axis_point = np.array(get_constriction(RCSB_ID))

    # Initialize the grid processor
    grid_processor = OptimizedVoxelGrid(use_gpu=True)

    # Get and transform points
    residues = filter_residues_parallel(ribosome_entities(RCSB_ID, 'R'), base_point, axis_point, R, H)
    points = np.array([atom.get_coord() for residue in residues for atom in residue.child_list])
    transformed_points = transform_points_optimized(points, base_point, axis_point)

    # Create the mask
    mask, (x, y, z) = grid_processor.create_point_cloud_mask(
        transformed_points,
        radius=R,
        height=H,
        voxel_size=Vsize,
        radius_around_point=ATOM_SIZE
    )

    # Process results
    points = np.where(~mask)
    empty_coordinates = np.column_stack((
        x[points[0]], 
        y[points[1]], 
        z[points[2]]
    ))
    world_coords = transform_points_from_C0_optimized(empty_coordinates, base_point, axis_point)
    
    # Visualization
    occupied_points = pv.PolyData(empty_coordinates)
    world_coords_pv = pv.PolyData(world_coords)
    return occupied_points, world_coords_pv

if __name__ == '__main__':
    main()