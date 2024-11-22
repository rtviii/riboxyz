import numpy as np
import pyvista as pv

def visualize_dot_projection(point, base_point, axis_point):
    """
    Visualizes how the dot product creates a projection onto the cylinder axis
    """
    # Create PyVista plotter
    plotter = pv.Plotter()
    
    # Calculate vectors
    point_vector = point - base_point
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length
    
    # Calculate projection
    projection = np.dot(point_vector, axis_unit)
    projection_point = base_point + projection * axis_unit
    
    # Add points
    plotter.add_mesh(pv.Sphere(radius=0.1, center=base_point), 
                    color='green', label='Base Point')
    plotter.add_mesh(pv.Sphere(radius=0.1, center=point), 
                    color='red', label='Test Point')
    plotter.add_mesh(pv.Sphere(radius=0.1, center=projection_point), 
                    color='blue', label='Projection Point')
    
    # Create and add vector arrows
    # Original vector (point_vector)
    arrow1 = pv.Arrow(start=base_point, direction=point_vector, 
                     tip_length=0.25, shaft_radius=0.05)
    plotter.add_mesh(arrow1, color='red', label='point_vector')
    
    # Axis vector
    arrow2 = pv.Arrow(start=base_point, direction=axis, 
                     tip_length=0.25, shaft_radius=0.05)
    plotter.add_mesh(arrow2, color='gray', label='cylinder axis')
    
    # Projection vector
    projection_vector = projection * axis_unit
    arrow3 = pv.Arrow(start=base_point, direction=projection_vector, 
                     tip_length=0.25, shaft_radius=0.05)
    plotter.add_mesh(arrow3, color='blue', label='projection')
    
    # Add right angle indicator if vectors aren't parallel
    if not np.allclose(point_vector, projection_vector):
        # Create right angle marker between projection and original point
        right_angle = create_right_angle_marker(projection_point, point)
        plotter.add_mesh(right_angle, color='yellow', label='Right Angle')
    
    # Add text showing the dot product value
    plotter.add_text(f"Dot product (projection length): {projection:.2f}",
                    position='upper_left', font_size=12)
    
    # Angle between vectors
    angle = np.arccos(np.clip(np.dot(point_vector/np.linalg.norm(point_vector), axis_unit), -1.0, 1.0))
    plotter.add_text(f"Angle between vectors: {np.degrees(angle):.2f}°",
                    position='upper_right', font_size=12)
    
    # Add coordinate system and legend
    plotter.show_axes()
    plotter.add_legend()
    plotter.view_isometric()
    plotter.show()

def create_right_angle_marker(proj_point, point, size=0.2):
    """Creates a right angle marker between the projection and the point"""
    # Create points for right angle symbol
    p1 = proj_point
    p2 = point
    direction = p2 - p1
    perp = np.cross(direction, [0, 0, 1])
    if np.allclose(perp, 0):
        perp = np.cross(direction, [0, 1, 0])
    perp = perp / np.linalg.norm(perp) * size
    
    points = np.array([
        p1,
        p1 + perp,
        p1 + perp + size * direction/np.linalg.norm(direction)
    ])
    
    lines = pv.PolyData()
    lines.points = points
    lines.lines = np.array([3, 0, 1, 2])
    
    return lines

# Example usage
base_point = np.array([0, 0, 0])
axis_point = np.array([0, 0, 5])
test_point = np.array([2, 2, 3])
visualize_dot_projection(test_point, base_point, axis_point)