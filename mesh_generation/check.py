
from mesh_generation.paths import alphashape_ensemble_LSU
import pyvista as pv

RCSB_ID="6Z6K"
plotter               = pv.Plotter()
mesh_  = pv.read(alphashape_ensemble_LSU(RCSB_ID))
plotter.add_mesh(mesh_, opacity=0.5)
plotter.show()