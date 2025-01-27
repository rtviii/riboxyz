from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from Bio.PDB.Entity import Entity
from Bio.PDB.MMCIFParser import FastMMCIFParser
import numpy as np
from Bio.PDB.Chain import Chain
import numpy as np
from typing import Callable, Literal, Tuple, List, TypeVar, Union, Optional
import numpy as np
import pyvista as pv
pv.global_theme.allow_empty_mesh = True


def create_cylinder_from_points(base_point, axis_point, radius, height):
    """
    Create a cylinder using two points: a base point and a point defining the axis direction.
    
    Parameters:
    base_point: array-like, [x, y, z] starting point of cylinder
    axis_point: array-like, [x, y, z] point defining cylinder axis direction
    radius: float, radius of cylinder
    height: float, total height of cylinder (can extend past axis_point)
    
    Returns:
    cylinder: PyVista cylinder mesh (triangulated)
    direction: numpy array of normalized direction vector
    """
    base_point = np.array(base_point)
    axis_point = np.array(axis_point)
    
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    ).triangulate()
    
    return cylinder, direction

def intersect_hull_with_cylinder(points, base_point, axis_point, radius, height):
    """
    Create convex hull from point cloud and intersect it with a cylinder 
    defined by two points.
    
    Parameters:
    points: numpy array of shape (N, 3) containing the point cloud
    base_point: array-like, [x, y, z] base point of cylinder
    axis_point: array-like, [x, y, z] point defining cylinder axis
    radius: float, radius of cylinder
    height: float, height of cylinder
    
    Returns:
    dict containing:
        - intersection: PyVista mesh of the intersection
        - hull: PyVista mesh of the convex hull
        - cylinder: PyVista mesh of the cylinder
        - direction: numpy array of cylinder axis direction
    """
    # Create convex hull
    cloud = pv.PolyData(points)
    hull = cloud.delaunay_3d().extract_surface()
    hull = hull.triangulate()
    
    # Create cylinder
    cylinder, direction = create_cylinder_from_points(
        base_point, axis_point, radius, height
    )
    
    # Compute intersection
    intersection = hull.boolean_intersection(cylinder)
    
    return {
        'intersection': intersection,
        'hull': hull,
        'cylinder': cylinder,
        'direction': direction
    }
    
def visualize_intersection(result, points, base_point, axis_point):
    """
    Visualize the intersection, hull, cylinder, and points.
    """
    plotter = pv.Plotter()
    
    # Add hull with transparency
    plotter.add_mesh(result['hull'], color='blue', opacity=0.3, label='Convex Hull')
    
    # Add cylinder with transparency
    plotter.add_mesh(result['cylinder'], color='red', opacity=0.3, label='Cylinder')
    
    # Add intersection
    plotter.add_mesh(result['intersection'], color='green', label='Intersection')
    
    # Add original points
    plotter.add_points(points, color='black', point_size=10, label='Points')
    
    # Add cylinder base and axis points
    plotter.add_points(np.array([base_point]), color='red', point_size=15, label='Cylinder Base')
    plotter.add_points(np.array([axis_point]), color='blue', point_size=15, label='Cylinder Axis Point')
    
    plotter.add_legend()
    return plotter

# -------------------------------------------
T = TypeVar('T')
