import os
import numpy as np
import pyvista as pv
from ribctl import EXIT_TUNNEL_WORK
from plotly import graph_objects as go

RCSB_ID   = "6Z6K"
file_path = os.path.join(EXIT_TUNNEL_WORK,"{}_poisson_recon.ply".format(RCSB_ID))  # Update with your file path

mesh    = pv.read(file_path)
plotter = pv.Plotter()
plotter.add_mesh(mesh, opacity=0.5)
# plotter.show()

# Add a sphere with specified size and color
points1 =np.array([[40,40,40], [50,50,50]])
# plotter.add_points(points1, render_points_as_spheres=True, point_size=100.0, color='red', style='points_gaussian')
plotter.add_points(points1,  point_size=100.0, color='red', style='points_gaussian')
points2 =np.array([[10,10,10], [100,100,100]])
plotter.add_points(points2, render_points_as_spheres=True, point_size=30.0, color='blue')

# Add a text label
plotter.add_text('Label Text', position='upper_left', font_size=18)

# Customize the view
plotter.show(auto_close=False)


# figb = go.Figure(data=[mesh])
# figb.update_layout(title_text="tit", title_x=0.5, width=800, height=800)
# figb.show()

