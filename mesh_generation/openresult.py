import pyvista as pv

RCSB_ID = "6Z6K"
# Load the PLY file
file_path = "{}_poisson_recon.ply".format(RCSB_ID)  # Update with your file path
mesh = pv.read(file_path)

# Plot the mesh
plotter = pv.Plotter()
plotter.add_mesh(mesh)
plotter.show()