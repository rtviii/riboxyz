import numpy as np
import pyvista as pv
from mesh_generation.mes_visualization import visualize_pointcloud
from mesh_generation.mesh_full_pipeline import expand_bbox_atoms_to_spheres

def create_cylinder_voxel_grid(
    base_point: np.ndarray, 
    axis_point: np.ndarray, 
    radius: float, 
    voxel_size: float
) -> tuple[np.ndarray, list]:
    """Creates a 3D binary voxel grid representing a cylinder"""
    # Calculate cylinder axis and height
    axis            = axis_point - base_point
    height          = np.linalg.norm(axis)
    axis_normalized = axis / height
    
    # Calculate grid dimensions for an axis-aligned bounding box
    padding = 2
    grid_dims = np.array([
        2 * int(np.ceil(radius / voxel_size)) + 2 * padding,
        2 * int(np.ceil(radius / voxel_size)) + 2 * padding,
        int(np.ceil(height / voxel_size)) + 2 * padding
    ])
    
    
    # Create grid
    voxel_grid = np.zeros(grid_dims, dtype=int)
    
    # Create a coordinate grid in the cylinder's local space
    z = np.arange(grid_dims[2]) * voxel_size
    y = (np.arange(grid_dims[1]) - grid_dims[1]//2) * voxel_size
    x = (np.arange(grid_dims[0]) - grid_dims[0]//2) * voxel_size
    
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    # Calculate which points are inside the cylinder
    radial_distance = np.sqrt(X**2 + Y**2)
    voxel_grid = (radial_distance <= radius) & (Z >= 0) & (Z <= height)
    
    # Store transforms
    transforms = [base_point, axis_normalized, voxel_size, grid_dims]
    
    return voxel_grid.astype(int), transforms

def world_to_grid_coords(
    points: np.ndarray, 
    base_point: np.ndarray, 
    axis_normalized: np.ndarray,
    grid_dims: np.ndarray,
    voxel_size: float
) -> np.ndarray:
    """Convert world coordinates to grid indices"""
    # Center points on base_point
    centered = points - base_point
    
    # Create rotation matrix to align cylinder with z-axis
    z_axis = np.array([0, 0, 1])
    if not np.allclose(axis_normalized, z_axis):
        v = np.cross(z_axis, axis_normalized)
        s = np.linalg.norm(v)
        c = np.dot(z_axis, axis_normalized)
        
        if s != 0:
            v = v / s
            V = np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])
            R = np.eye(3) + s * V + (1 - c) * (V @ V)
            centered = centered @ R.T
    
    # Convert to grid coordinates
    grid_center = np.array([grid_dims[0]//2, grid_dims[1]//2, 0])
    indices = np.floor(centered / voxel_size + grid_center).astype(int)
    
    # Filter valid indices
    valid = np.all((indices >= 0) & (indices < grid_dims), axis=1)
    return indices[valid]

def visualize_cylinder_voxel_grid(voxel_grid: np.ndarray, transforms: list, show_edges: bool = True):
    """Visualize the voxel grid"""
    base_point, axis_normalized, voxel_size, _ = transforms
    
    # Create grid
    grid = pv.ImageData()
    grid.dimensions = np.array(voxel_grid.shape) + 1
    grid.spacing = (voxel_size, voxel_size, voxel_size)
    grid.cell_data["values"] = voxel_grid.flatten(order="F")
    
    # Convert to unstructured grid for transformation
    grid = grid.cast_to_unstructured_grid()
    threshold = grid.threshold(0.5)
    outline = grid.outline()
    
    # Create rotation matrix
    z_axis = np.array([0, 0, 1])
    if not np.allclose(axis_normalized, z_axis):
        v = np.cross(z_axis, axis_normalized)
        s = np.linalg.norm(v)
        c = np.dot(z_axis, axis_normalized)
        
        if s != 0:
            v = v / s
            V = np.array([[0, -v[2], v[1]],
                         [v[2], 0, -v[0]],
                         [-v[1], v[0], 0]])
            R = np.eye(3) + s * V + (1 - c) * (V @ V)
            
            # Create 4x4 transformation matrix
            transform = np.eye(4)
            transform[:3, :3] = R
            transform[:3, 3] = base_point
            
            threshold.transform(transform)
            outline.transform(transform)
    else:
        threshold.translate(base_point)
        outline.translate(base_point)
    
    # Create plotter
    plotter = pv.Plotter()
    plotter.set_background("white")
    
    # Add meshes
    plotter.add_mesh(
        threshold,
        color="lightblue",
        opacity=0.6,
        show_edges=show_edges,
        edge_color="darkblue",
        line_width=1,
    )
    plotter.add_mesh(outline, color="black", line_width=2, opacity=1)
    plotter.add_axes()
    
    plotter.show()

if __name__ == "__main__":
    # Your input parameters
    base_point = np.array([179.15499878, 179.46508789, 160.99293518])
    axis_point = np.array([199.4345495, 264.11536385, 51.34317183])
    axis_point = np.array([230.4345495, 294.11536385, 81.34317183])
    radius     = 40
    voxel_size = 1

    # 1. Create the cylinder grid
    voxel_grid, transforms = create_cylinder_voxel_grid(base_point, axis_point, radius, voxel_size)
    
    # 2. Load your point cloud
    bbox_atoms_expanded = np.load("bbox_atoms_expanded.npy")
    visualize_pointcloud(bbox_atoms_expanded)
    exit()
    
    # 3. Invert the grid
    voxel_grid = 1 - voxel_grid
    
    # 4. Convert points to grid coordinates and set them
    valid_indices = world_to_grid_coords(
        bbox_atoms_expanded,
        transforms[0],  # base_point
        transforms[1],  # axis_normalized
        transforms[3],  # grid_dims
        transforms[2]   # voxel_size
    )
    
    # 5. Set points in grid
    if len(valid_indices) > 0:
        voxel_grid[valid_indices[:, 0], valid_indices[:, 1], valid_indices[:, 2]] = 1
    
    # 6. Visualize
    visualize_cylinder_voxel_grid(voxel_grid, transforms, show_edges=True)
