import argparse
from pprint import pprint
import subprocess
import typing
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import json
import os
import numpy as np
import numpy as np
from mesh_generation.mesh_bbox_extraction import (
    encode_atoms,
    open_tunnel_csv,
    parse_struct_via_bbox,
    parse_struct_via_centerline,
)
from compas.geometry import bounding_box
from mesh_generation.mesh_paths import *
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
import warnings
warnings.filterwarnings("ignore")


def apply_poisson_reconstruction(surf_estimated_ptcloud_path: str, output_path: str, recon_depth:int=6, recon_pt_weight:int=3):
    import plyfile
    # The documentation can be found at https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version16.04/ in "PoissonRecon" binary
    command = [
        POISSON_RECON_BIN,
        "--in",
        surf_estimated_ptcloud_path,
        "--out",
        output_path,
        "--depth",
        str(recon_depth),
        "--pointWeight",
        str(recon_pt_weight),
        "--threads 8"
    ]
    process = subprocess.run(command, capture_output=True, text=True)

    if process.returncode == 0:
        print(">>PoissonRecon executed successfully.")
        print(">>Wrote {}".format(output_path))
        # Convert the plyfile to asciii

        data = plyfile.PlyData.read(output_path)
        data.text = True
        ascii_duplicate =output_path.split(".")[0] + "_ascii.ply"
        data.write(ascii_duplicate)
        print(">>Wrote {}".format(ascii_duplicate))
    else:
        print(">>Error:", process.stderr)

def ptcloud_convex_hull_points(pointcloud: np.ndarray, ALPHA:float, TOLERANCE:float) -> np.ndarray:
    assert pointcloud is not None
    cloud       = pv.PolyData(pointcloud)
    grid        = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points

def ptcloud_convex_hull_points_and_gif(pointcloud: np.ndarray, ALPHA: float, TOLERANCE: float, output_path: str, n_frames: int = 180) -> np.ndarray:

    from tqdm import tqdm
    from PIL import Image
    assert pointcloud is not None
    
    # Create the convex hull
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface()
    
    # Set up the plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(convex_hull, color='lightblue', show_edges=True)
    plotter.add_points(cloud, color='red', point_size=5)
    
    # Set up the camera
    plotter.camera_position = 'xy'
    
    # Create frames
    frames = []
    for i in tqdm(range(n_frames)):
        plotter.camera.azimuth = i * (360 / n_frames)
        plotter.render()
        image = plotter.screenshot(transparent_background=False, return_img=True)
        frames.append(Image.fromarray(image))
    
    # Save as GIF
    frames[0].save(output_path, save_all=True, append_images=frames[1:], duration=50, loop=0)
    
    print(f"GIF saved to {output_path}")
    
    return convex_hull.points

def estimate_normals(convex_hull_surface_pts: np.ndarray, output_path: str, kdtree_radius=None, kdtree_max_nn=None, correction_tangent_planes_n=None): 
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=kdtree_radius, max_nn=kdtree_max_nn) )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    o3d.io.write_point_cloud(output_path, pcd)
    return pcd

from PIL import Image
import os

def estimate_normals_and_create_gif(convex_hull_surface_pts: np.ndarray, output_path: str, kdtree_radius=None, kdtree_max_nn=None, correction_tangent_planes_n=None):
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=kdtree_radius, max_nn=kdtree_max_nn))
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
        
    # Set up the visualizer
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.add_geometry(pcd)
    opt = vis.get_render_option()
    opt.point_size = 5
    opt.point_show_normal = True
    
    # Set up the camera
    ctr = vis.get_view_control()
    ctr.set_zoom(0.8)
    
    # Capture frames
    frames = []
    for i in range(360):
        ctr.rotate(10.0, 0.0)  # Rotate 10 degrees around the z-axis
        vis.update_geometry(pcd)
        vis.poll_events()
        vis.update_renderer()
        image = vis.capture_screen_float_buffer(False)
        frames.append(Image.fromarray((np.array(image) * 255).astype(np.uint8)))
    
    vis.destroy_window()
    
    frames[0].save(output_path, save_all=True, append_images=frames[1:], duration=50, loop=0)
    
    print(f"GIF saved to {output_path}")