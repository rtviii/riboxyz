import numpy as np
import pyvista as pv
force_float=False

def explore_3d_grid():
    # Create non-symmetrical coordinate arrays
    x = np.linspace(0,11,10)        # 3 points, unevenly spaced
    y = np.linspace(0,11,10)     # 4 points
    z = np.linspace(0,11,10)        # 3 points, big spacing
    
    print("Coordinate arrays:")
    print("x:", x, "  (size:", len(x), ")")
    print("y:", y, "  (size:", len(y), ")")
    print("z:", z, "  (size:", len(z), ")")
    
    # Create 3D grid
    X, Y, Z = np.meshgrid(x, y, z, indexing='ij')
    
    print("\nGrid shapes:")
    print("X shape:", X.shape)  # Should be (3, 4, 3)
    
    # Create a sample voxel grid (just for demonstration)
    # Let's make a simple pattern where voxel is True if x+y+z > 5
    voxel_grid = (X + Y + Z) > 5
    
    # Let's look at some specific voxels
    print("\nExploring specific voxels:")
    
    # Let's look at a few example points
    example_indices = [
        (0, 0, 0),  # First corner
        (2, 3, 2),  # Last corner
        (1, 2, 1),  # Some middle point
        (2, 1, 1),  # Some middle point

    ]
    
    for i, j, k in example_indices:
        print(f"\nVoxel at index ({i}, {j}, {k}):")
        print(f"  x coordinate: {x[i]}")
        print(f"  y coordinate: {y[j]}")
        print(f"  z coordinate: {z[k]}")
        print(f"  Position in 3D space: ({x[i]}, {y[j]}, {z[k]})")
        print(f"  Voxel value: {voxel_grid[i,j,k]}")
    
    # Visualize occupied voxels
    occupied = np.where(voxel_grid)
    points = np.column_stack((
        x[occupied[0]], 
        y[occupied[1]], 
        z[occupied[2]]
    ))
    
    # Also create a point cloud of ALL grid points to show the structure
    all_points = np.column_stack((
        X.ravel(), Y.ravel(), Z.ravel()
    ))
    
    plotter = pv.Plotter()
    
    # Add all grid points (small, gray)
    point_cloud_all = pv.PolyData(all_points)
    plotter.add_mesh(point_cloud_all, color='gray', point_size=5, label='All Grid Points')
    
    # Add occupied points (larger, blue)
    point_cloud_occupied = pv.PolyData(points)
    plotter.add_mesh(point_cloud_occupied, color='blue', point_size=15, label='Occupied Points')
    
    # Add some coordinate labels for selected points
    for i, j, k in example_indices:
        pos = np.array([x[i], y[j], z[k]])
        plotter.add_point_labels(
            [pos], 
            [f"({i},{j},{k})\n({x[i]:.1f}, {y[j]:.1f}, {z[k]:.1f})"],
            point_size=20,
        )
    
    plotter.add_legend()
    plotter.show_axes()
    plotter.show()

explore_3d_grid()