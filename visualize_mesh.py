
import sys
import pyvista as pv
import os

def visualize_mesh(mesh_id):
    mesh_path = os.path.join('wenjun_data_meshes', f'{mesh_id}_NPET_MESH.ply')
    
    if not os.path.exists(mesh_path):
        print(f"Mesh file not found: {mesh_path}")
        return
    
    # Create a new plotter
    plotter = pv.Plotter()
    
    # Load and add the mesh
    mesh = pv.read(mesh_path)
    plotter.add_mesh(mesh, color='lightblue', show_edges=True)
    
    # Set a title with the ID
    plotter.add_text(f"Ribosome ID: {mesh_id}", position='upper_left', font_size=12)
    
    # Show the plotter in a separate window
    plotter.show()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        visualize_mesh(sys.argv[1])
    else:
        print("Please provide a mesh ID as an argument")
