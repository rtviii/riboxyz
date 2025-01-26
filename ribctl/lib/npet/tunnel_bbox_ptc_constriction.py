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


def ribosome_entities(rcsb_id:str, cifpath:str, level=Literal['R']|Literal[ 'A' ])->list[Entity]:
    structure = FastMMCIFParser(QUIET=True).get_structure(rcsb_id, cifpath)
    residues = []
    [residues.extend(chain) for chain in structure.child_list[0] ]
    if level == 'R':
        return residues
    elif level == 'A':
        atoms = []
        [atoms.extend(a) for a in residues]
        return atoms
    else:
        raise

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

def make_cylinder_predicate(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray] 
) -> Callable[[T], bool]:
    def predicate(obj: T) -> bool:
        position = position_getter(obj)
        return is_point_in_cylinder(position, base_point, axis_point, radius, height)
    return predicate

def get_residue_position(residue):
    return residue.center_of_mass()
# def get_residue_position(residue):
#     return residue.get_coord()

def _worker_process_chunk(chunk_data):
    """
    Worker function that processes a chunk of residues.
    Takes a tuple of (residue_positions, base_point, axis_point, radius, height, indices)
    Returns indices of residues that are inside the cylinder.
    """
    positions, base_point, axis_point, radius, height, indices = chunk_data
    
    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(pos, base_point, axis_point, radius, height):
            results.append(indices[i])
    
    return results

def filter_residues_parallel(
    residues: List[T],
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None
) -> List[T]:
    """
    Filter residues in parallel using ProcessPoolExecutor.
    
    Parameters:
    -----------
    residues : List[T]
        List of residue objects to filter
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    chunk_size : Optional[int]
        Size of chunks to process in parallel. If None, calculated automatically
    max_workers : Optional[int]
        Maximum number of worker processes. If None, uses CPU count
    
    Returns:
    --------
    List[T]
        Filtered list of residues whose positions lie within the cylinder
    """
    # Set defaults for parallel processing parameters
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    if chunk_size is None:
        # Aim for ~4 chunks per worker for better load balancing
        chunk_size = max(1, len(residues) // (max_workers * 4))
    
    # Pre-compute all positions and create index mapping
    positions = np.array([get_residue_position(r) for r in residues])
    
    indices = list(range(len(residues)))
    index_chunks = [
        indices[i:i + chunk_size] 
        for i in range(0, len(indices), chunk_size)
    ]
    
    # Create data chunks for processing
    chunks_data = [ (positions[idx], base_point, axis_point, radius, height, idx) for idx in index_chunks ]
    
    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))
    
    # Flatten results and get corresponding residues
    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]

def visualize_filtered_residues(
    filtered_residues: List[T],
    all_residues: Optional[List[T]],  # If None, won't show unfiltered residues
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray] = lambda x: x.center_of_mass(),
    point_size: float = 5,
    opacity: float = 0.3,
    show: bool = True,
    screenshot_path: Optional[str] = None,
    window_size: tuple = (1024, 768)
) -> Optional[pv.Plotter]:
    """
    Visualize filtered residues alongside the cylinder that was used for filtering.
    
    Parameters:
    -----------
    filtered_residues : List[T]
        List of residues that passed the cylinder filter
    all_residues : Optional[List[T]]
        Complete list of residues before filtering. If provided, will show
        unfiltered residues in gray
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    position_getter : Callable[[T], np.ndarray]
        Function to extract position from residue object
    point_size : float
        Size of points in visualization
    opacity : float
        Opacity of cylinder (0.0-1.0)
    show : bool
        Whether to display the plot immediately
    screenshot_path : Optional[str]
        If provided, saves screenshot to this path
    window_size : tuple
        Size of the visualization window (width, height)
        
    Returns:
    --------
    Optional[pv.Plotter]
        Returns plotter if show=False, None otherwise
    """
    # Initialize plotter
    plotter = pv.Plotter(window_size=window_size)
    
    # Get positions of filtered residues
    filtered_positions = np.array([
        position_getter(res) for res in filtered_residues
    ])
    
    # Create and add cylinder
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    )
    
    # Add cylinder with wireframe and solid style
    plotter.add_mesh(
        cylinder, 
        style='wireframe', 
        color='red', 
        line_width=2, 
        label='Cylinder'
    )
    plotter.add_mesh(
        cylinder, 
        style='surface', 
        color='red', 
        opacity=opacity,
    )
    
    # If all_residues provided, show unfiltered residues in gray
    if all_residues is not None:
        all_positions = np.array([
            position_getter(res) for res in all_residues
        ])
        plotter.add_points(
            all_positions, 
            color='gray', 
            point_size=point_size-2, 
            opacity=0.3,
            label='Unfiltered Residues'
        )
    
    # Add filtered points
    plotter.add_points(
        filtered_positions, 
        color='blue', 
        point_size=point_size,
        label='Filtered Residues'
    )
    
    # Add base and axis points for reference
    plotter.add_points(
        np.array([base_point]), 
        color='green', 
        point_size=point_size*2, 
        label='Base Point'
    )
    plotter.add_points(
        np.array([axis_point]), 
        color='yellow', 
        point_size=point_size*2, 
        label='Axis Point'
    )
    
    # Add axis line
    axis_line = pv.Line(base_point, axis_point)
    plotter.add_mesh(
        axis_line, 
        color='white', 
        line_width=2, 
        label='Cylinder Axis'
    )
    
    # Add legend
    plotter.add_legend()
    
    # Set camera position for better initial view
    plotter.camera_position = 'iso'
    plotter.enable_eye_dome_lighting()  # Improves point visibility
    
    # Add text with stats
    stats_text = (
        f'Total filtered residues: {len(filtered_residues)}\n'
        f'Cylinder height: {height:.1f}\n'
        f'Cylinder radius: {radius:.1f}'
    )
    if all_residues:
        stats_text = f'Total residues: {len(all_residues)}\n' + stats_text
        
    plotter.add_text(
        stats_text,
        position='upper_left',
        font_size=12,
        color='white'
    )
        
    if show:
        plotter.show()
        return None
        
    return plotter

def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float

) -> bool:

    point      = np.asarray(point)
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)
    
    # Calculate cylinder axis direction vector
    axis        = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit   = axis / axis_length
    
    # Calculate vector from base to point
    point_vector = point - base_point
    
    # Project vector onto cylinder axis
    projection = np.dot(point_vector, axis_unit)
    
    # Calculate perpendicular vector from axis to point
    projection_point = base_point + projection * axis_unit
    radial_vector    = point - projection_point
    
    # Calculate radial distance
    radial_distance = np.linalg.norm(radial_vector)
    
    # Check if point is inside cylinder
    return (radial_distance <= radius) and (0 <= projection <= height)

# def get_npet_cylinder_residues(rcsb_id:str, radius, height):
#     residues   = ribosome_entities(rcsb_id, 'R')
#     base_point = np.array(PTC_location(rcsb_id).location)
#     axis_point = np.array( get_constriction(rcsb_id) )

#     return filter_residues_parallel(
#         residues,
#         base_point,
#         axis_point,
#         radius,
#         height,
#     ),  base_point, axis_point


