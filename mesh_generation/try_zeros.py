import numpy as np
import pyvista as pv


coords = np.array([[0,0,0],[1,1,1],[2,2,2]])
shape = (5,5,5)
grid = np.zeros(shape)
grid[
    coords[:,0],
    coords[:,1],
    coords[:,2]
] = 1


indices = np.asanyarray(np.where(grid==1))

pl                        = pv.Plotter()
_                         = pl.add_points(grid[:,:-1,2:], color='blue', point_size=10, 
                                        #   render_points_as_spheres=True,
                                          show_edges=True)
_ = pl.add_axes(line_width=5,     
    cone_radius=0.6,
    shaft_length=0.7,
    tip_length=0.3,
    ambient=0.5,
    label_size=(0.4, 0.16),
    )
pl.show()
# _                         = pl.add_points(pts, color='r', point_size=4)
# _                         = pl.add_points(main_cluster, opacity=0.5, color='b' ,point_size=2)
# pl.add_text('ALPHA VAL: {}'.format(8), position='upper_left', font_size=20, shadow=True, font='courier', color='black')