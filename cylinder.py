import os
from typing import Tuple
import numpy as np

from mesh_generation.mes_visualization import visualize_pointcloud
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import get_npet_cylinder_residues


import numpy as np
import pyvista as pv
from typing import Tuple


def visualize_voxel_grid(
    voxel_grid: np.ndarray,
    coordinates: Tuple[np.ndarray, np.ndarray, np.ndarray],
    transformed_base,
    transformed_axis,
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
    
    plotter.add_mesh(pv.Sphere(radius=4, center=transformed_axis), color='blue', label='Axis Point')
    plotter.add_mesh(pv.Sphere(radius=4, center=transformed_base), color='red', label='Base Point')
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
points = np.array([atom.get_coord() for residue in residues for atom in residue.child_list])
voxel_size = 1

translation, rotation = get_transformation_to_C0(base, axis)
t_base = ( base + translation ) @ rotation.T
t_axis = ( axis + translation ) @ rotation.T

if os.path.exists('points.npy'):
    points = np.load('points.npy')
else:
    points = np.array([ residue.center_of_mass() for residue  in residues])
    np.save('points.npy', points)
    print("SAVED")


nx = ny = int(2 * radius / voxel_size) + 1
nz = int(height / voxel_size) + 1
x = np.linspace(-radius, radius, nx)
y = np.linspace(-radius, radius, ny)
z = np.linspace(0, height, nz)
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')


transformed = transform_points_to_C0(points, base, axis)
X_I = np.round(transformed[:,0])
Y_I = np.round(transformed[:,1])
Z_I = np.round(transformed[:,2])

cylinder_mask = (np.sqrt(X**2 + Y**2) <= radius)
hollow_cylinder = ~cylinder_mask

# 3. Create point cloud mask
point_cloud_mask = np.zeros_like(X, dtype=bool)
for point in zip(X_I, Y_I, Z_I):
    point_cloud_mask |= (np.abs(X - point[0])<=2) & (np.abs(Y == point[1])<2) & (np.abs(Z == point[2])<2)

final_mask = hollow_cylinder | point_cloud_mask

occupied = np.where(final_mask)
points = np.column_stack((
    x[occupied[0]], 
    y[occupied[1]], 
    z[occupied[2]]
))
occupied_points = pv.PolyData(points)
visualize_pointcloud(occupied_points)