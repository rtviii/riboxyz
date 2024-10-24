from ribctl.etl.etl_assets_ops import RibosomeOps
import numpy as np
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from scipy.spatial.distance import cdist

from ribctl.lib.libtax import Taxid

def midpoint(p1, p2):
    return (p1 + p2) / 2

def find_closest_points(points1, points2):
    """
    Find the pair of points (one from each set) that are closest to each other.
    
    Parameters:
    points1: numpy array of shape (N, 3) for first set of 3D points
    points2: numpy array of shape (M, 3) for second set of 3D points
    
    Returns:
    tuple: (point1, point2, distance) where point1 and point2 are the closest points
           and distance is their Euclidean distance
    """
    # Convert inputs to numpy arrays if they aren't already
    points1 = np.asarray(points1)
    points2 = np.asarray(points2)
    
    # Calculate pairwise distances between all points
    distances = cdist(points1, points2)
    
    # Find the indices of the minimum distance
    i, j = np.unravel_index(distances.argmin(), distances.shape)
    
    # Get the closest points and their distance
    closest_point1 = points1[i]
    closest_point2 = points2[j]
    min_distance = distances[i, j]
    
    return closest_point1, closest_point2

def bbox_axis(ptc:np.ndarray, constriction:np.ndarray):
    ...

def rrna_ptcloud(rcsb_id:str)->np.ndarray:
    ro = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile().mitochondrial
    if is_mitochondrial:
        lsurna = ro.get_poly_by_polyclass('mt16SrRNA')
    else:
        lsurna = ro.get_poly_by_polyclass('23SrRNA')
        if lsurna is None:
            lsurna = ro.get_poly_by_polyclass('25SrRNA')
        elif lsurna is None:
            lsurna = ro.get_poly_by_polyclass('28SrRNA')
    if lsurna is None:
        raise ValueError("Could not find 23SrRNA, 25SrRNA, 28SrRNA, or mt16SrRNA in {}".format(rcsb_id))

    structure = ro.biopython_structure()
    rna_chain:Chain = structure[0][lsurna.auth_asym_id]
    return np.array([r.center_of_mass()  for r in rna_chain.child_list])


def get_constriction(rcsb_id: str):
    ro = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile().mitochondrial
    if is_mitochondrial:
        uL4   = ro.get_poly_by_polyclass('uL4m')
        uL22 = ro.get_poly_by_polyclass('uL22m')
    else:
        uL4   = ro.get_poly_by_polyclass('uL4')
        uL22 = ro.get_poly_by_polyclass('uL22')

    if uL4 is None or uL22 is None:
        raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))

    structure = ro.biopython_structure()

    uL4_c       :Chain = structure[0][uL4.auth_asym_id]
    uL22_c      :Chain = structure[0][uL22.auth_asym_id]

    uL4_coords  = [(r.center_of_mass() ) for r in uL4_c.child_list]
    uL22_coords = [(r_.center_of_mass() ) for r_ in uL22_c.child_list]

    return midpoint(*find_closest_points(uL4_coords, uL22_coords))

   


    
    



# get_constriction('3J9M')
# rrna_ptcloud('3J9M')

import numpy as np
import pyvista as pv


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

# Example usage
# if __name__ == "__main__":
#     # Parameters you'll need to provide:
#     # 1. Point cloud
#     points = rrna_ptcloud('3J9M')
    
#     # 2. Cylinder parameters
#     base_point = np.array([0., 0., 0.])  # Replace with your base point
#     axis_point = np.array([1., 1., 1.])  # Replace with your axis point
#     radius = 1.0                         # Replace with your radius
#     height = 3.0                         # Replace with your height
    
#     try:
#         # Compute intersection
#         result = intersect_hull_with_cylinder(
#             points=points,
#             base_point=base_point,
#             axis_point=axis_point,
#             radius=radius,
#             height=height
#         )
        
#         # Print some information about the results
#         print("Intersection Properties:")
#         print(f"Volume: {result['intersection'].volume}")
#         print(f"Surface Area: {result['intersection'].area}")
#         print(f"Number of points: {len(result['intersection'].points)}")
        
#         # Visualize
#         plotter = visualize_intersection(result, points, base_point, axis_point)
#         plotter.show()
        
#     except Exception as e:
#         print(f"Error occurred: {str(e)}")
#         print("\nDebug information:")
#         print(f"Number of points: {len(points)}")
#         print(f"Distance between base and axis points: "
#               f"{np.linalg.norm(axis_point - base_point)}")
