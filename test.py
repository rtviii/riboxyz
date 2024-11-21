import numpy as np
import pyvista as pv

from ribctl.lib.npet.tunnel_bbox_ptc_constriction import get_npet_cylinder_residues


def create_cylinder_voxel_grid(
    base_point: np.ndarray, axis_point: np.ndarray, radius: float, voxel_size: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Creates a 3D binary voxel grid representing a cylinder.
    The grid is always aligned with the cylinder's axis.
    """
    # Calculate cylinder properties
    axis = axis_point - base_point
    height = np.linalg.norm(axis)
    axis_normalized = axis / height

    # Calculate grid dimensions to fully contain cylinder
    grid_radius = int(np.ceil(radius / voxel_size))
    grid_height = int(np.ceil(height / voxel_size))

    # Create grid with padding
    pad = 2
    grid_dims = np.array(
        [2 * grid_radius + 2 * pad, 2 * grid_radius + 2 * pad, grid_height + 2 * pad]
    )

    # Create coordinate grid in cylinder's local space
    x = (np.arange(-grid_radius - pad, grid_radius + pad)) * voxel_size
    y = (np.arange(-grid_radius - pad, grid_radius + pad)) * voxel_size
    z = (np.arange(-pad, grid_height + pad)) * voxel_size

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")

    # Calculate radial distance in local space (cylinder is vertical here)
    radial_distance = np.sqrt(X**2 + Y**2)

    # Create cylinder mask
    voxel_grid = (radial_distance <= radius) & (Z >= 0) & (Z <= height)

    # Store transformation info
    transform_vectors = [base_point, axis_normalized, voxel_size]

    return voxel_grid.astype(int), transform_vectors


def visualize_cylinder_voxel_grid(
    voxel_grid: np.ndarray, transform_vectors: np.ndarray, show_edges: bool = True
) -> None:
    """
    Visualize the voxel grid with the cylinder inside using PyVista.
    Handles proper transformation to world space.
    """
    base_point, axis_normalized, voxel_size = transform_vectors

    # Create initial grid
    grid = pv.ImageData()
    grid.dimensions = np.array(voxel_grid.shape) + 1
    grid.spacing = (voxel_size, voxel_size, voxel_size)
    grid.cell_data["values"] = voxel_grid.flatten(order="F")

    # Convert to UnstructuredGrid for transformation
    grid = grid.cast_to_unstructured_grid()
    threshold = grid.threshold(0.5)
    outline = grid.outline()

    # Create the rotation matrix to transform from local to world space
    z_axis = np.array([0, 0, 1])
    if not np.allclose(axis_normalized, z_axis):
        # Calculate rotation matrix to align z-axis with cylinder axis
        v = np.cross(z_axis, axis_normalized)
        s = np.linalg.norm(v)
        c = np.dot(z_axis, axis_normalized)

        if s != 0:
            v = v / s
            V = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            R = np.eye(3) + s * V + (1 - c) * (V @ V)

            # Create 4x4 transformation matrix
            transform = np.eye(4)
            transform[:3, :3] = R
            transform[:3, 3] = base_point

            # Apply transformation
            threshold.transform(transform)
            outline.transform(transform)
    else:
        # Just translate if axis is already vertical
        threshold.translate(base_point)
        outline.translate(base_point)

    # Create plotter
    plotter = pv.Plotter()
    plotter.set_background("white")

    # Add cylinder voxels
    plotter.add_mesh(
        threshold,
        color="lightblue",
        opacity=0.6,
        show_edges=show_edges,
        edge_color="darkblue",
        line_width=1,
    )

    # Add transformed bounding box
    plotter.add_mesh(outline, color="black", line_width=2, opacity=1)

    # Add axes for reference
    plotter.add_axes()

    plotter.show()


def calculate_axis_endpoint(
    base_point: np.ndarray, point_on_axis: np.ndarray, height: float
) -> np.ndarray:
    """
    Calculate the endpoint of a cylinder axis given base point, any point on desired axis, and height.

    Parameters:
    -----------
    base_point : np.ndarray
        Starting point of cylinder (x,y,z)
    point_on_axis : np.ndarray
        Any point that lies on the desired axis line
    height : float
        Desired height of the cylinder

    Returns:
    --------
    axis_point : np.ndarray
        The endpoint of the cylinder axis
    """
    # Calculate direction vector
    direction = point_on_axis - base_point

    # Normalize direction vector
    direction_normalized = direction / np.linalg.norm(direction)

    # Calculate axis endpoint by adding height-sized vector in direction
    axis_point = base_point + direction_normalized * height

    return axis_point


if __name__ == "__main__":

    RCSB_ID = "5AFI"
    npet_residues, bp, ap, radius, height = get_npet_cylinder_residues(RCSB_ID)
    atoms_poss = np.array(
        [atom.get_coord() for residue in npet_residues for atom in residue.child_list]
    )

    # Test with vertical cylinder
    # base_point = np.array([0.0, 0.0, 0.0])
    # axis_point = np.array([0.0, 0.0, 40.0])

    # base_point =bp
    # axis_point =calculate_axis_endpoint (bp,ap,140)
    base_point = np.array([179.15499878, 179.46508789, 160.99293518])
    axis_point = np.array([199.4345495, 264.11536385, 51.34317183])

    print("Base point: ", base_point)
    print("axis point: ", axis_point)
    radius = 40
    voxel_size = 1

    grid, transforms = create_cylinder_voxel_grid(
        base_point, axis_point, radius, voxel_size
    )
    visualize_cylinder_voxel_grid(grid, transforms, show_edges=True)

    # # Test with tilted cylinder
    # base_point = np.array([0.0, 0.0, 0.0])
    # axis_point = np.array([20.0, 20.0, 20.0])

    # grid, transforms = create_cylinder_voxel_grid(
    #     base_point, axis_point, radius, voxel_size
    # )
    # visualize_cylinder_voxel_grid(grid, transforms, show_edges=True)
