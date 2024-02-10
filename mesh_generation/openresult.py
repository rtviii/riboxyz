import pyvista as pv

# read the data
grid = pv.read('6Z6K_tunnel_delaunay.vtk')

# plot the data with an automatically created Plotter
grid.plot(show_scalar_bar=False, show_axes=False)