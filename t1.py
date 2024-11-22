import os
import numpy as np
from typing import Tuple, List, Optional
import pyvista as pv
import numpy as np
from typing import Tuple, List, Optional
import pyvista as pv

from mesh_generation.mes_visualization import visualize_pointcloud
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import get_npet_cylinder_residues

def visualize_pyvista_split(
        points         : np.ndarray,
        cylinder_points: np.ndarray,
        empty_voxels   : np.ndarray,
        base_point     : np.ndarray,
        axis_point     : np.ndarray,
        radius         : float,
        height         : float,
        grid_params    : dict) -> None: 
    """
    Visualize results using PyVista with split views.
    
    Args:

        points         : original point cloud
        cylinder_points: points inside cylinder
        empty_voxels   : coordinates of empty voxels
        base_point     : cylinder base center
        axis_point     : point to determine cylinder axis direction
        radius         : cylinder radius
        height         : cylinder height
        grid_params    : grid parameters from create_cylinder_voxel_grid
    """
    pl = pv.Plotter(shape=(1, 2))
    
    # First subplot: Original points and cylinder
    pl.subplot(0, 0)
    pl.add_title("Original Points and Cylinder")
    
    # Add original points
    point_cloud = pv.PolyData(points)
    pl.add_mesh(point_cloud, color='blue', point_size=5, 
               render_points_as_spheres=True, label='Original Points')

    
    # Add cylinder points
    if len(cylinder_points) > 0:
        cylinder_cloud = pv.PolyData(cylinder_points)
        pl.add_mesh(cylinder_cloud, color='red', point_size=5,
                   render_points_as_spheres=True, label='Points in Cylinder')
    
    # Add cylinder outline
    axis = grid_params['axis']
    top_point = grid_params['top_point']
    cylinder = pv.Cylinder(center=(base_point + top_point)/2,
                         direction=axis,
                         radius=radius,
                         height=height)

    pl.add_mesh(cylinder, style='wireframe', color='black', 
               opacity=0.5, label='Cylinder')
    
    pl.add_axes()
    pl.add_legend()
    
    # Second subplot: Voxel grid and empty spaces
    pl.subplot(0, 1)
    pl.add_title("Cylinder Voxel Grid and Empty Space")
    
    # Add cylinder outline
    pl.add_mesh(cylinder, style='wireframe', color='black', 
               opacity=0.5, label='Cylinder')
    
    # Add empty voxels
    if len(empty_voxels) > 0:
        empty_cloud = pv.PolyData(empty_voxels)
        pl.add_mesh(empty_cloud, color='green', point_size=5,
                    opacity=0.5, 
                   label='Empty Voxels')
    
    # Add cylinder points for reference
    if len(cylinder_points) > 0:
        pl.add_mesh(cylinder_cloud, color='red', point_size=5,
                    label='Points in Cylinder')
    
    pl.add_axes()
    pl.add_legend()
    
    # Link cameras
    pl.link_views()
    
    # Show the plot
    pl.show()

def compute_cylinder_params(base_point: np.ndarray, 
                          axis_point: np.ndarray, 
                          radius: float, 
                          height: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute cylinder parameters including axis direction and rotation matrix.
    """
    # Compute axis direction
    axis = axis_point - base_point
    axis = axis / np.linalg.norm(axis)
    
    # Compute actual top point using height
    top_point = base_point + axis * height
    
    # Create rotation matrix to align cylinder with z-axis
    z_axis        = np.array([0, 0, 1])
    rotation_axis = np.cross(axis, z_axis)

    if np.all(rotation_axis == 0):
        rotation_matrix = np.eye(3)
    else:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        angle = np.arccos(np.dot(axis, z_axis))
        
        # Rodriguez rotation formula
        K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                     [rotation_axis[2], 0, -rotation_axis[0]],
                     [-rotation_axis[1], rotation_axis[0], 0]])
        rotation_matrix = (np.eye(3) + np.sin(angle) * K + 
                         (1 - np.cos(angle)) * np.matmul(K, K))
    
    return axis, top_point, rotation_matrix

def transform_to_cylinder_coords(points: np.ndarray,
                               base_point: np.ndarray,
                               rotation_matrix: np.ndarray) -> np.ndarray:
    """
    Transform points from original to cylinder coordinates.
    """
    return np.dot(points - base_point, rotation_matrix)

def transform_from_cylinder_coords(points: np.ndarray,
                                 base_point: np.ndarray,
                                 rotation_matrix: np.ndarray) -> np.ndarray:
    """
    Transform points from cylinder coordinates back to original coordinates.
    """
    return np.dot(points, rotation_matrix.T) + base_point

def points_in_cylinder(points: np.ndarray, 
                      base_point: np.ndarray, 
                      axis_point: np.ndarray, 
                      radius: float,
                      height: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract points that lie within the cylinder.
    Returns both mask and transformed points.
    """
    axis, top_point, rotation_matrix = compute_cylinder_params(
        base_point, axis_point, radius, height)
    
    # Transform points to cylinder coordinate system
    transformed_points = transform_to_cylinder_coords(points, base_point, rotation_matrix)
    
    # Check radial distance
    xy_dist = np.sqrt(transformed_points[:, 0]**2 + transformed_points[:, 1]**2)
    radial_mask = xy_dist <= radius
    
    # Check height
    height_mask = (transformed_points[:, 2] >= 0) & (transformed_points[:, 2] <= height)
    
    mask = radial_mask & height_mask
    return mask, transformed_points[mask]

def create_cylinder_voxel_grid(transformed_points: np.ndarray, radius: float, height: float, voxel_size: float) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    # Calculate grid dimensions
    nx = ny = int(2 * radius / voxel_size)  # Remove +1 as we'll center the grid
    nz = int(height / voxel_size)
    
    # Create voxel coordinates (cell centers)
    # Center the grid at (0,0) in xy-plane, start at 0 in z
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)
    
    # Initialize grid
    voxel_grid = np.zeros((nx, ny, nz), dtype=bool)
    
    # Find occupied voxels
    ix = np.clip(np.searchsorted(x, transformed_points[:, 0]) - 1, 0, nx-1)
    iy = np.clip(np.searchsorted(y, transformed_points[:, 1]) - 1, 0, ny-1)
    iz = np.clip(np.searchsorted(z, transformed_points[:, 2]) - 1, 0, nz-1)
    
    voxel_grid[ix, iy, iz] = True
    return voxel_grid, (x, y, z)

def get_empty_voxel_coordinates(voxel_grid: np.ndarray, 
                              grid_coords: Tuple[np.ndarray, np.ndarray, np.ndarray],
                              radius: float,
                              base_point: np.ndarray,
                              rotation_matrix: np.ndarray) -> np.ndarray:
    """
    Extract coordinates of empty voxels in original coordinate system.
    """
    x, y, z = grid_coords
    empty_voxels = []
    
    for i in range(len(x)):
        for j in range(len(y)):
            for k in range(len(z)):
                if not voxel_grid[i, j, k]:
                    point = np.array([x[i], y[j], z[k]])
                    if np.sqrt(point[0]**2 + point[1]**2) <= radius:
                        original_point = transform_from_cylinder_coords( point, base_point, rotation_matrix)
                        empty_voxels.append(original_point)
    
    return np.array(empty_voxels) if empty_voxels else np.array([])

def main():
    # Load data
    RCSB_ID = '3J7Z'
    base_point = np.array([179.15499878, 179.46508789, 160.99293518])
    axis_point = np.array([199.4345495, 264.11536385, 51.34317183])
    radius = 40
    height = 150
    voxel_size = 1





    # points = np.load("bbox_atoms_expanded.npy")
    if os.path.exists('points.npy'):
        points = np.load('points.npy')
    else:
        npet_residues, bp,ap,radius, height       = get_npet_cylinder_residues(RCSB_ID, radius, height)
        points = np.array([r.center_of_mass() for r in npet_residues])

    visualize_pointcloud(points)
    # axis, top_point, rotation_matrix = compute_cylinder_params( base_point, axis_point, radius, height)
    
    # Find points in cylinder
    # mask, transformed_points = points_in_cylinder(
    #     points, base_point, axis_point, radius, height)
    # cylinder_points = points[mask]
    
    # # Create voxel grid in cylinder coordinates
    # voxel_grid, grid_coords = create_cylinder_voxel_grid( transformed_points, radius, height, voxel_size)
    
    # empty_voxels = get_empty_voxel_coordinates( voxel_grid, grid_coords, radius, base_point, rotation_matrix)
    
    # # Visualize results
    # visualize_pyvista_split(points, cylinder_points, empty_voxels, base_point, axis_point, radius, height, {'axis': axis, 'top_point': top_point})

if __name__ == "__main__":
    main()