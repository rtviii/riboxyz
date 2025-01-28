import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
import numpy as np
from scipy.spatial import ConvexHull
import alphashape
import trimesh
from typing import List, Tuple
import multiprocessing as mp
from functools import partial

from ribctl.asset_manager.asset_types import AssetType

warnings.filterwarnings("ignore")

def cif_to_point_cloud(cif_path: str) -> np.ndarray:
    """Optimized point cloud generation using numpy operations"""
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    # Pre-allocate approximate size based on typical protein structure
    coords = np.zeros((50000, 3))  # Adjust size based on your typical structures
    idx = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if idx >= coords.shape[0]:
                        coords = np.vstack((coords, np.zeros((10000, 3))))
                    coords[idx] = atom.coord
                    idx += 1
    
    return coords[:idx]

def batch_sample_points(args: Tuple[trimesh.Trimesh, np.ndarray, np.ndarray, int]) -> np.ndarray:
    """Worker function for parallel point sampling"""
    mesh, bbox_min, bbox_max, batch_size = args
    points = np.random.uniform(low=bbox_min, high=bbox_max, size=(batch_size, 3))
    mask = mesh.contains(points)
    return points[mask]

def sample_within_alpha_shape_parallel(alpha_shape_mesh: trimesh.Trimesh, 
                                     num_samples: int,
                                     batch_size: int = 10000,
                                     num_processes: int = None) -> np.ndarray:
    """
    Optimized parallel sampling within alpha shape using batch processing
    """
    if num_processes is None:
        num_processes = max(1, mp.cpu_count() - 1)
    
    bbox_min, bbox_max = alpha_shape_mesh.bounds
    
    # Estimate required points based on volume ratio
    hull = ConvexHull(alpha_shape_mesh.vertices)
    volume_ratio = alpha_shape_mesh.volume / hull.volume
    estimated_points = int(num_samples / volume_ratio * 1.2)  # Add 20% buffer
    
    # Prepare batches for parallel processing
    num_batches = max(1, estimated_points // batch_size)
    args = [(alpha_shape_mesh, bbox_min, bbox_max, batch_size)] * num_batches
    
    # Use parallel processing to generate points
    sampled_points = []
    with mp.Pool(num_processes) as pool:
        for points in pool.imap_unordered(batch_sample_points, args):
            sampled_points.extend(points)
            if len(sampled_points) >= num_samples:
                break
    
    # Ensure exact number of points
    sampled_points = np.array(sampled_points[:num_samples])
    return sampled_points

def produce_alpha_contour(RCSB_ID: str, alpha: float, num_samples: int = 5000) -> None:
    """
    Optimized alpha shape generation with progress tracking
    """
    print(f"Processing {RCSB_ID} with alpha={alpha}")
    
    # Generate point cloud
    cifpath = AssetType.MMCIF.get_path(RCSB_ID)
    point_cloud = cif_to_point_cloud(cifpath)
    print(f"Generated point cloud with {len(point_cloud)} points")
    
    # Create initial alpha shape
    print("Constructing initial alpha shape...")
    alpha_shape = alphashape.alphashape(point_cloud, alpha)
    
    # Get largest component
    components = alpha_shape.split(only_watertight=False)
    alpha_shape_largest = max(components, key=lambda c: abs(c.volume))
    print(f"Largest component volume: {alpha_shape_largest.volume:.2f}")
    
    # Parallel resampling
    print(f"Resampling {num_samples} points...")
    new_points = sample_within_alpha_shape_parallel(
        alpha_shape_largest,
        num_samples,
        batch_size=10000
    )
    
    # Generate final alpha shape
    print("Constructing final alpha shape...")
    alpha_shape_renew = alphashape.alphashape(new_points, alpha)
    
    # Save result
    path = AssetType.ALPHA_SHAPE.get_path(RCSB_ID)
    alpha_shape_renew.export(path, file_type="ply", encoding="ascii")
    print(f"Saved alpha shape to {path}")