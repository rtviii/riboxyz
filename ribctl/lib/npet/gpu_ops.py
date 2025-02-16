import cupy as cp
import cupyx.scipy.spatial as cu_spatial
from cupyx.scipy.spatial.distance import cdist as gpu_cdist
import numpy as np
import open3d as o3d
import torch
from torch_geometric.nn import knn_graph
import pytorch3d.ops as ops

def gpu_point_cloud_processing(coordinates: np.ndarray) -> np.ndarray:
    """
    GPU-accelerated version of coordinate processing
    """
    # Transfer data to GPU
    coords_gpu = cp.asarray(coordinates)
    
    # Parallel center of mass calculation
    if coords_gpu.ndim == 3:  # If processing residues with multiple atoms
        centers = cp.mean(coords_gpu, axis=1)
    else:
        centers = coords_gpu
        
    return cp.asnumpy(centers)

def gpu_normal_estimation(points: np.ndarray, k: int = 30) -> o3d.geometry.PointCloud:
    """
    GPU-accelerated normal estimation using PyTorch
    """
    # Convert to torch tensor
    device = torch.device('cuda')
    points_torch = torch.from_numpy(points).float().to(device)
    
    # Compute KNN graph
    edge_index = knn_graph(points_torch, k=k)
    
    # Get neighboring points for each point
    neighbors = points_torch[edge_index[1]].view(-1, k, 3)
    points_expanded = points_torch.unsqueeze(1).expand(-1, k, -1)
    
    # Compute covariance matrices
    centered = neighbors - points_expanded
    covariance = torch.bmm(centered.transpose(1, 2), centered)
    
    # Compute eigenvectors (normals are smallest eigenvector)
    eigenvalues, eigenvectors = torch.linalg.eigh(covariance)
    normals = eigenvectors[:, :, 0]  # Smallest eigenvector
    
    # Convert back to numpy and create Open3D point cloud
    points_np = points_torch.cpu().numpy()
    normals_np = normals.cpu().numpy()
    
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(points_np)
    pcd.normals = o3d.utility.Vector3dVector(normals_np)
    
    return pcd

def gpu_dbscan(points: np.ndarray, eps: float, min_samples: int) -> np.ndarray:
    """
    GPU-accelerated DBSCAN implementation
    """
    # Transfer data to GPU
    points_gpu = cp.asarray(points)
    
    # Compute pairwise distances on GPU
    distances = gpu_cdist(points_gpu, points_gpu)
    
    # Find neighbors
    neighbors = cp.where(distances <= eps, 1, 0)
    neighbor_counts = cp.sum(neighbors, axis=1)
    
    # Core points
    core_points = cp.where(neighbor_counts >= min_samples)[0]
    
    # Cluster assignment (simplified implementation)
    labels = cp.full(len(points), -1)
    current_cluster = 0
    
    for point_idx in core_points:
        if labels[point_idx] != -1:
            continue
            
        # Assign new cluster
        cluster_points = cp.where(neighbors[point_idx] == 1)[0]
        labels[cluster_points] = current_cluster
        
        # Expand cluster
        to_check = cluster_points
        while len(to_check) > 0:
            new_points = cp.unique(cp.concatenate([
                cp.where(neighbors[p] == 1)[0] for p in to_check
            ]))
            new_points = new_points[labels[new_points] == -1]
            labels[new_points] = current_cluster
            to_check = new_points
            
        current_cluster += 1
    
    return cp.asnumpy(labels)

def gpu_surface_reconstruction(points: np.ndarray, normals: np.ndarray) -> torch.Tensor:
    """
    GPU-accelerated surface reconstruction using PyTorch3D
    """
    device = torch.device('cuda')
    points_torch = torch.from_numpy(points).float().to(device)
    normals_torch = torch.from_numpy(normals).float().to(device)
    
    # Use PyTorch3D's built-in reconstruction
    mesh = ops.sample_points_from_meshes(
        points_torch.unsqueeze(0),
        normals_torch.unsqueeze(0)
    )
    
    return mesh

# Example usage in alpha_contour_via_poisson_recon:
def gpu_accelerated_alpha_contour(rcsb_id: str, verbose: bool = False):
    # ... (previous setup code remains the same)
    
    # GPU-accelerated point cloud processing
    ptcloud = gpu_point_cloud_processing(cif_to_point_cloud(...))
    
    # GPU-accelerated surface points calculation
    surface_pts = quick_surface_points(ptcloud, d3d_alpha, d3d_tol, d3d_offset)
    
    # GPU-accelerated normal estimation
    normal_estimated_pcd = gpu_normal_estimation(
        surface_pts,
        k=max_nn
    )
    
    # GPU-accelerated surface reconstruction
    mesh = gpu_surface_reconstruction(
        np.asarray(normal_estimated_pcd.points),
        np.asarray(normal_estimated_pcd.normals)
    )
    
    return mesh