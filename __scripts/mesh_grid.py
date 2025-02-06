import pyvista as pv
import os
from ribctl.asset_manager.asset_types import AssetType

def visualize_multiple_meshes(mesh_paths: list, rcsb_ids: list|None=None, gif: bool=False, gif_name: str|None=None):
    # Validate inputs
    if len(mesh_paths) > 24:
        raise ValueError("Maximum number of meshes is 24")
    
    if rcsb_ids and len(rcsb_ids) != len(mesh_paths):
        raise ValueError("Number of RCSB IDs must match number of mesh paths")
        
    shape = (4, 6)
    plotter = pv.Plotter(shape=shape, off_screen=gif)
    
    n_meshes = len(mesh_paths)
    
    for i in range(n_meshes):
        row = i // 6
        col = i % 6
        
        plotter.subplot(row, col)
        
        if os.path.exists(mesh_paths[i]):
            _ = plotter.add_mesh(pv.read(mesh_paths[i]), opacity=0.8)
        else:
            box = pv.Box(bounds=(-1, 1, -1, 1, -1, 1))
            _ = plotter.add_mesh(box, style='wireframe', opacity=0.3, color='gray')
            
        plotter.add_axes(
            line_width=1,
            cone_radius=0.5,
            shaft_length=0.5,
            tip_length=0.2,
            ambient=0.5,
            label_size=(0.15, 0.6)
        )
        
        if rcsb_ids:
            plotter.add_text(
                'RCSB_ID:{}'.format(rcsb_ids[i] if rcsb_ids else ""),
                position='upper_right',
                font_size=10,
                shadow=True,
                font='courier',
                color='black'
            )
            
        plotter.show_grid(
            n_xlabels=4,
            n_ylabels=4,
            n_zlabels=4,
            font_size=6
        )

    if gif:
        plotter.open_gif(gif_name)
        for angle in range(0, 360, 2):
            for i in range(n_meshes):
                row = i // 6
                col = i % 6
                plotter.subplot(row, col)
                plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {gif_name}")
    else:
        plotter.show()

structs = ["5X8T", "3J9M", "5NJT", "5AFI", "5JVG", "5MYJ", "5O60", "5V7Q", "5NRG", 
           "4Y4P", "4V9F", "4V6U", "6EK0", "5T2A", "3J79", "5GAK", "4V7E", "5T5H", 
           "5XXB", "5XY3", "3J7Z", "5VP2", "4UG0"]
meshpaths = list(map(lambda x: AssetType.NPET_MESH.get_path(x), structs))
visualize_multiple_meshes(meshpaths, structs)
