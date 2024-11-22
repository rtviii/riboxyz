import os
from typing import Tuple
import numpy as np

from mesh_generation.mes_visualization import visualize_pointcloud
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import get_npet_cylinder_residues


import numpy as np
import pyvista as pv
from typing import Tuple

def create_cylinder_voxel_grid(
    radius: float,
    height: float,
    voxel_size: float
) -> Tuple[np.ndarray, Tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """
    Creates a voxel grid representation of a cylinder centered at origin, aligned with z-axis.
    
    Args:
        radius: Radius of cylinder
        height: Height of cylinder
        voxel_size: Size of each voxel
    
    Returns:
        Tuple containing:
        - 3D boolean array representing the voxel grid
        - Tuple of coordinate arrays (x, y, z) for the voxel centers
    """
    # Calculate grid dimensions
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    
    # Create grid coordinates centered at origin for x,y
    # and starting at 0 for z
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)
    
    # Create 3D coordinate grid
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Calculate radial distance for each point
    R = np.sqrt(X**2 + Y**2)
    
    # Create voxel grid
    voxel_grid = (R <= radius)
    
    return voxel_grid, (x, y, z)

def visualize_voxel_grid(
    voxel_grid: np.ndarray,
    coordinates: Tuple[np.ndarray, np.ndarray, np.ndarray],
    filled_points: np.ndarray = np.ndarray([])
):
    plotter = pv.Plotter()
    x, y, z = coordinates
    
    occupied = np.where(voxel_grid)
    points   = np.column_stack((
        x[occupied[0]], 
        y[occupied[1]], 
        z[occupied[2]]
    ))
    
    # Create point cloud
    point_cloud = pv.PolyData(points)
    plotter.add_mesh(point_cloud, color='black', point_size=2, opacity=0.1,show_edges=True,edge_color='gray')
    plotter.add_mesh(filled_points, color='red', point_size=6, opacity=0.5, show_edges=True,edge_color='gray')
    
    plotter.show_axes()
    plotter.show_grid()
    plotter.show()

import numpy as np
import pyvista as pv
from typing import Tuple

def get_transformation_to_C0(base_point: np.ndarray, axis_point: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute transformation matrices to move arbitrary cylinder to C0 configuration.
    Returns translation vector and rotation matrix.
    """
    # Get cylinder axis vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length
    
    # Get rotation that aligns axis_unit with [0, 0, 1]
    z_axis = np.array([0, 0, 1])
    
    # Use Rodrigues rotation formula to find rotation matrix
    # that rotates axis_unit to z_axis
    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])  # 180-degree rotation around x-axis
    else:
        v = np.cross(axis_unit, z_axis)
        s = np.linalg.norm(v)
        c = np.dot(axis_unit, z_axis)
        v_skew = np.array([[0, -v[2], v[1]],
                          [v[2], 0, -v[0]],
                          [-v[1], v[0], 0]])
        R = np.eye(3) + v_skew + (v_skew @ v_skew) * (1 - c) / (s * s)
    
    return -base_point, R

def transform_points_to_C0(points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray) -> np.ndarray:
    """
    Transform points from arbitrary cylinder configuration to C0.
    """
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    
    points_translated  = points + translation
    points_transformed = points_translated @ rotation.T
    
    return points_transformed


RCSB_ID = '3J7Z'
radius     = 40
height     = 80
residues, base, axis = get_npet_cylinder_residues(RCSB_ID, radius=radius, height=height)

points = np.array([a.get_coord() for residue in residues for a in residue])
atom_points =  points
np.save('points.npy', atom_points)

# base_point = np.array([179.15499878, 179.46508789, 160.99293518])
# axis_point = np.array([199.4345495, 264.11536385, 51.34317183])
voxel_size = 1

# # points = np.load("bbox_atoms_expanded.npy")
# if os.path.exists('points.npy'):
#     points = np.load('points.npy')
# else:
#     npet_residues, bp,ap,radius, height       = get_npet_cylinder_residues(RCSB_ID, radius, height)
#     points = np.array([r.center_of_mass() for r in npet_residues])

grid,(x,y,z) = create_cylinder_voxel_grid(radius, height,1)
visualize_voxel_grid(grid,(x,y,z,), points)
transformed= transform_points_to_C0(points, base, axis)
visualize_voxel_grid(grid,(x,y,z,), transformed)


