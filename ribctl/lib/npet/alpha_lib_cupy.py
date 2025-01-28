import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
import numpy as np
import cupy as cp
import alphashape
import trimesh
from typing import Union, Tuple
from time import time

from ribctl.asset_manager.asset_types import AssetType

warnings.filterwarnings("ignore")

def cif_to_point_cloud(cif_path: str) -> np.ndarray:
    """Optimized point cloud generation with pre-allocation"""
    parser = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    coords = np.zeros((50000, 3), dtype=np.float32)  # Use float32 for GPU compatibility
    idx = 0
    
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if idx >= coords.shape[0]:
                        coords = np.vstack((coords, np.zeros((10000, 3), dtype=np.float32)))
                    coords[idx] = atom.coord
                    idx += 1
    
    return coords[:idx]

def points_in_mesh_gpu(points: cp.ndarray, 
                      vertices: cp.ndarray, 
                      faces: cp.ndarray) -> cp.ndarray:
    """
    CUDA kernel for checking if points are inside the mesh using ray casting
    """
    def get_ray_direction() -> cp.ndarray:
        """Generate a random ray direction that's not aligned with any face normal"""
        while True:
            direction = cp.random.randn(3).astype(cp.float32)
            direction /= cp.linalg.norm(direction)
            normals = cp.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]], vertices[faces[:, 2]] - vertices[faces[:, 0]])
            normals /= cp.linalg.norm(normals, axis=1, keepdims=True)
            if not cp.any(cp.abs(cp.dot(normals, direction)) < 1e-6):
                return direction
    
    ray_direction = get_ray_direction()
    
    def ray_triangle_intersect(origins: cp.ndarray) -> cp.ndarray:
        """Vectorized ray-triangle intersection test"""
        v0 = vertices[faces[:, 0]]
        v1 = vertices[faces[:, 1]]
        v2 = vertices[faces[:, 2]]
        
        # Compute triangle normals and areas
        normals = cp.cross(v1 - v0, v2 - v0)
        areas = cp.linalg.norm(normals, axis=1) / 2
        normals /= 2 * areas[:, None]
        
        # Expand origins for broadcasting
        origins = origins[:, None, :]
        
        # Ray-plane intersection
        d = cp.sum(normals * (v0 - origins), axis=2)
        t = d / cp.dot(normals, ray_direction)
        
        # Intersection points
        p = origins + t[:, :, None] * ray_direction
        
        # Barycentric coordinates
        edge0 = v1 - v0
        edge1 = v2 - v0
        v2_v0 = v2 - v0
        p_v0 = p - v0
        
        dot00 = cp.sum(edge0 * edge0, axis=1)
        dot01 = cp.sum(edge0 * edge1, axis=1)
        dot02 = cp.sum(edge0 * p_v0, axis=2)
        dot11 = cp.sum(edge1 * edge1, axis=1)
        dot12 = cp.sum(edge1 * p_v0, axis=2)
        
        invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01)
        u = (dot11 * dot02 - dot01 * dot12) * invDenom
        v = (dot00 * dot12 - dot01 * dot02) * invDenom
        
        # Check if intersection point is inside triangle
        mask = (u >= 0) & (v >= 0) & (u + v <= 1) & (t > 0)
        return cp.sum(mask, axis=1) % 2 == 1

    return ray_triangle_intersect(points)

def sample_within_alpha_shape_gpu(alpha_shape_mesh: trimesh.Trimesh,
                                num_samples: int,
                                batch_size: int = 100000) -> np.ndarray:
    """
    GPU-accelerated point sampling within alpha shape
    """
    # Transfer mesh data to GPU
    vertices = cp.array(alpha_shape_mesh.vertices, dtype=cp.float32)
    faces = cp.array(alpha_shape_mesh.faces, dtype=cp.int32)
    
    # Get bounding box
    bbox_min = cp.min(vertices, axis=0)
    bbox_max = cp.max(vertices, axis=0)
    
    sampled_points = []
    start_time = time()
    
    while len(sampled_points) < num_samples:
        # Generate random points on GPU
        points = cp.random.uniform(
            bbox_min, bbox_max,
            size=(batch_size, 3)
        ).astype(cp.float32)
        
        # Check which points are inside the mesh
        inside_mask = points_in_mesh_gpu(points, vertices, faces)
        
        # Get points inside the mesh
        inside_points = cp.asnumpy(points[inside_mask])
        sampled_points.extend(inside_points)
        
        if len(sampled_points) >= num_samples:
            print(f"Sampled {len(sampled_points)} points in {time() - start_time:.2f} seconds")
            break
    
    return np.array(sampled_points[:num_samples])

def produce_alpha_contour_cupy(RCSB_ID: str, alpha: float, num_samples: int = 5000) -> None:
    """
    GPU-accelerated alpha shape generation
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
    
    # GPU-accelerated resampling
    print(f"Resampling {num_samples} points using GPU...")
    new_points = sample_within_alpha_shape_gpu(
        alpha_shape_largest,
        num_samples,
        batch_size=100000
    )
    
    # Generate final alpha shape
    print("Constructing final alpha shape...")
    alpha_shape_renew = alphashape.alphashape(new_points, alpha)
    
    # Save result
    path = AssetType.ALPHA_SHAPE.get_path(RCSB_ID)
    alpha_shape_renew.export(path, file_type="ply", encoding="ascii")
    print(f"Saved alpha shape to {path}")

# # Example usage
# if __name__ == "__main__":
#     produce_alpha_contour_cupy("1J5E", alpha=12.0, num_samples=5000)