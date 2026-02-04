from pprint import pprint
from typing import Callable, List, Optional
from matplotlib import pyplot as plt
import pyvista as pv
import json
import numpy as np
import open3d as o3d

from ribctl.lib.npet.kdtree_approach import T

hexcolors = {
    'aliceblue'           : '#F0F8FF',
    'antiquewhite'        : '#FAEBD7',
    'aquamarine'          : '#7FFFD4',
    'azure'               : '#F0FFFF',
    'beige'               : '#F5F5DC',
    'bisque'              : '#FFE4C4',
    'black'               : '#000000',
    'blanchedalmond'      : '#FFEBCD',
    'blue'                : '#0000FF',
    'blueviolet'          : '#8A2BE2',
    'brown'               : '#A52A2A',
    'burlywood'           : '#DEB887',
    'cadetblue'           : '#5F9EA0',
    'chartreuse'          : '#7FFF00',
    'chocolate'           : '#D2691E',
    'coral'               : '#FF7F50',
    'cornflowerblue'      : '#6495ED',
    'cornsilk'            : '#FFF8DC',
    'crimson'             : '#DC143C',
    'cyan'                : '#00FFFF',
    'darkblue'            : '#00008B',
    'darkcyan'            : '#008B8B',
    'darkgoldenrod'       : '#B8860B',
    'darkgray'            : '#A9A9A9',
    'darkgreen'           : '#006400',
    'darkkhaki'           : '#BDB76B',
    'darkmagenta'         : '#8B008B',
    'darkolivegreen'      : '#556B2F',
    'darkorange'          : '#FF8C00',
    'darkorchid'          : '#9932CC',
    'darkred'             : '#8B0000',
    'darksalmon'          : '#E9967A',
    'darkseagreen'        : '#8FBC8F',
    'darkslateblue'       : '#483D8B',
    'darkslategray'       : '#2F4F4F',
    'darkturquoise'       : '#00CED1',
    'darkviolet'          : '#9400D3',
    'deeppink'            : '#FF1493',
    'deepskyblue'         : '#00BFFF',
    'dimgray'             : '#696969',
    'dodgerblue'          : '#1E90FF',
    'firebrick'           : '#B22222',
    'floralwhite'         : '#FFFAF0',
    'forestgreen'         : '#228B22',
    'gainsboro'           : '#DCDCDC',
    'ghostwhite'          : '#F8F8FF',
    'gold'                : '#FFD700',
    'goldenrod'           : '#DAA520',
    'gray'                : '#808080',
    'green'               : '#008000',
    'greenyellow'         : '#ADFF2F',
    'honeydew'            : '#F0FFF0',
    'hotpink'             : '#FF69B4',
    'indianred'           : '#CD5C5C',
    'indigo'              : '#4B0082',
    'ivory'               : '#FFFFF0',
    'khaki'               : '#F0E68C',
    'lavender'            : '#E6E6FA',
    'lavenderblush'       : '#FFF0F5',
    'lawngreen'           : '#7CFC00',
    'lemonchiffon'        : '#FFFACD',
    'lightblue'           : '#ADD8E6',
    'lightcoral'          : '#F08080',
    'lightcyan'           : '#E0FFFF',
    'lightgoldenrodyellow': '#FAFAD2',
    'lightgray'           : '#D3D3D3',
    'lightgreen'          : '#90EE90',
    'lightpink'           : '#FFB6C1',
    'lightsalmon'         : '#FFA07A',
    'lightseagreen'       : '#20B2AA',
    'lightskyblue'        : '#87CEFA',
    'lightslategray'      : '#778899',
    'lightsteelblue'      : '#B0C4DE',
    'lightyellow'         : '#FFFFE0',
    'lime'                : '#00FF00',
    'limegreen'           : '#32CD32',
    'linen'               : '#FAF0E6',
    'magenta'             : '#FF00FF',
    'maroon'              : '#800000',
    'mediumaquamarine'    : '#66CDAA',
    'mediumblue'          : '#0000CD',
    'mediumorchid'        : '#BA55D3',
    'mediumpurple'        : '#9370DB',
    'mediumseagreen'      : '#3CB371',
    'mediumslateblue'     : '#7B68EE',
    'mediumspringgreen'   : '#00FA9A',
    'mediumturquoise'     : '#48D1CC',
    'mediumvioletred'     : '#C71585',
    'midnightblue'        : '#191970',
    'mintcream'           : '#F5FFFA',
    'mistyrose'           : '#FFE4E1',
    'moccasin'            : '#FFE4B5',
    'navajowhite'         : '#FFDEAD',
    'navy'                : '#000080',
    'oldlace'             : '#FDF5E6',
    'olive'               : '#808000',
    'olivedrab'           : '#6B8E23',
    'orange'              : '#FFA500',
    'orangered'           : '#FF4500',
    'orchid'              : '#DA70D6',
    'palegoldenrod'       : '#EEE8AA',
    'palegreen'           : '#98FB98',
    'paleturquoise'       : '#AFEEEE',
    'palevioletred'       : '#DB7093',
    'papayawhip'          : '#FFEFD5',
    'paraview_background' : '#52576e',
    'peachpuff'           : '#FFDAB9',
    'peru'                : '#CD853F',
    'pink'                : '#FFC0CB',
    'plum'                : '#DDA0DD',
    'powderblue'          : '#B0E0E6',
    'purple'              : '#800080',
    'raw_sienna'          : '#965434',
    'rebeccapurple'       : '#663399',
    'red'                 : '#FF0000',
    'rosybrown'           : '#BC8F8F',
    'royalblue'           : '#4169E1',
    'saddlebrown'         : '#8B4513',
    'salmon'              : '#FA8072',
    'sandybrown'          : '#F4A460',
    'seagreen'            : '#2E8B57',
    'seashell'            : '#FFF5EE',
    'sienna'              : '#A0522D',
    'silver'              : '#C0C0C0',
    'skyblue'             : '#87CEEB',
    'slateblue'           : '#6A5ACD',
    'slategray'           : '#708090',
    'snow'                : '#FFFAFA',
    'springgreen'         : '#00FF7F',
    'steelblue'           : '#4682B4',
    'tan'                 : '#D2B48C',
    'teal'                : '#008080',
    'thistle'             : '#D8BFD8',
    'tomato'              : '#FF6347',
    'turquoise'           : '#40E0D0',
    'violet'              : '#EE82EE',
    'wheat'               : '#F5DEB3',
    'white'               : '#FFFFFF',
    'whitesmoke'          : '#F5F5F5',
    'yellow'              : '#FFFF00',
    'yellowgreen'         : '#9ACD32',
    'tab:blue'            : '#1f77b4',
    'tab:orange'          : '#ff7f0e',
    'tab:green'           : '#2ca02c',
    'tab:red'             : '#d62728',
    'tab:purple'          : '#9467bd',
    'tab:brown'           : '#8c564b',
    'tab:pink'            : '#e377c2',
    'tab:gray'            : '#7f7f7f',
    'tab:olive'           : '#bcbd22',
    'tab:cyan'            : '#17becf',
}

diagram_tunnels = {
    "bacteria": [
        "4W29", # GOOD
        "6WD4", # GOOD
        "7UNV",
        "5DM6",

        "5O60",
        "8BUU",
        "6HMA",
        "7RYH", #GOOD

        "7MSZ",
        "7P7T",
        "5MYJ", # GOOD, weird structure
        "7JIL",
        # 3j9y -- good
    ],
    "eukaryota": [
        "6P5N",
        "7QGG",
        "4UG0", # GOOD
        "4U3M",

        "7OYB",
        "7OLC",
        "8EUI",
        "5XXB",

        "4V91",
        "6XU8",
        "4V7E",
        "8P5D",

        "8BTR",
        "3JBO",
        "7CPU",
        "7Q08",
        "6AZ3", #GOOD

        "5T5H",
        "5XY3",
        "7QEP",
    ],
    "archaea": ["4V6U",  # GOOD
                "4V9F",] # GOOD
}

FONT                  = 'courier'
CHAIN_PT_SIZE         = 8
PTC_PT_SIZE           = 40
CHAIN_LANDMARK_COLORS = ["purple","orange", "cornflowerblue", "cornsilk", "crimson", "darkblue", "darkcyan", "darkgoldenrod", "darkgray", "darkgreen", "darkkhaki", "darkmagenta", "darkolivegreen", "darkorange", "darkorchid", "darkred",
"rebeccapurple",
"rosybrown",
"royalblue",
"saddlebrown",
"salmon",
"sandybrown",
"seagreen"]
dbscan_pairs = [
    (2,33),
    (2.3, 57),
    (3, 123),
    (3.2, 145),
    (3.5, 175),
    (4.2, 280),
    (5, 490),
    (5.5,600)
]


def visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
    dbscan_cluster_dict: dict[int, list], eps, min_nbrs, base_point, axis_point, 
    refined, radius, height, output_path: str | None = None):
    """
    Visualize DBSCAN clustering results and optionally save to file.
    
    Parameters:
        dbscan_cluster_dict: Dictionary of clusters
        eps: Epsilon parameter used for DBSCAN
        min_nbrs: Minimum samples parameter used for DBSCAN
        base_point: Base point of cylinder
        axis_point: Axis point of cylinder
        refined: Refined points
        radius: Cylinder radius
        height: Cylinder height
        output_path: Path to save the visualization image
    """
    plotter = pv.Plotter(off_screen=output_path is not None)
    plotter.subplot(0, 0)

    # Print cluster sizes
    cluster_sizes = {k: len(v) for k, v in dbscan_cluster_dict.items()}
    for k, v in cluster_sizes.items():
        print(f"Cluster {k} has {v} points.")

    # Find the largest non-noise cluster
    largest_cluster = max(
        ((k, len(v)) for k, v in dbscan_cluster_dict.items() if k != -1),
        key=lambda x: x[1]
    )[0]

    # Generate color palette for regular clusters
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]

    # Create cylinder
    direction = np.array(axis_point) - np.array(base_point)
    direction = direction / np.linalg.norm(direction)
    center = np.array(base_point) + (direction * height/2)
    
    cylinder = pv.Cylinder(
        center=center,
        direction=direction,
        radius=radius,
        height=height
    )
    
    plotter.add_mesh(cylinder, opacity=0.1, color='gray', label='Cylinder', style='wireframe')

    # Add each cluster separately
    for cluster_label, coordinates in dbscan_cluster_dict.items():
        points = np.array(coordinates)
        if cluster_label == -1:
            # Noise points
            plotter.add_points(
                points,
                color='gray',
                opacity=0.1,
                point_size=1,
                label='Noise',
                style='points_gaussian',
                emissive=True
            )
        elif cluster_label == largest_cluster:
            # Largest cluster
            plotter.add_points(
                points,
                color='cyan',
                point_size=5,
                opacity=0.15,
                label=f'Cluster {cluster_label} (Largest)',
                style='points',
                render_points_as_spheres=True
            )
        else:
            # Regular clusters
            color = clusters_palette[(cluster_label * 2) % len(clusters_palette)]
            plotter.add_points(
                points,
                color='gray',
                opacity=0.3,
                point_size=2,
                label=f'Cluster {cluster_label}',
                style='points'
            )

    plotter.add_points(
        refined,
        color='blue',
        point_size=3,
        opacity=1,
        label='refined_points',
        style='points',
        render_points_as_spheres=True
    )
    
    # Add reference points
    plotter.add_points(
        np.array([base_point]),
        color='red',
        point_size=20,
        label='Base Point',
        style='points',
        render_points_as_spheres=True
    )
    plotter.add_points(
        np.array([axis_point]),
        color='red',
        point_size=20,
        label='Axis Point',
        style='points',
        render_points_as_spheres=True
    )

    # Add text information
    plotter.add_text(
        f'eps: {eps}\nmin_nbrs: {min_nbrs}',
        position='upper_left',
        font_size=20,
        shadow=True,
        font=FONT,
        color='black'
    )

    if output_path:
        plotter.screenshot(output_path)
        print(f"Saved visualization to {output_path}")
    else:
        plotter.show()
    
    return plotter

def DBSCAN_CLUSTERS_visualize_largest(positive_space: np.ndarray, dbscan_cluster_dict: dict[int, list], selected_cluster: np.ndarray, gif:bool=False, gif_name:str|None=None):
    plotter               = pv.Plotter(shape=(1, 2), off_screen=True)
    plotter.subplot(0,0)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 8)
    #? Visualize all clusters
    for k, v in dbscan_cluster_dict.items():
        print("Cluster {} has {} points.".format(k, len(v)))
    clusters_palette = dict(zip(range(-1, 60), plt.cm.terrain(np.linspace(0, 1, 60))))
    for k, v in clusters_palette.items():
        clusters_palette[k] = [*v[:3], 0.5]
    combined_cluster_colors = []
    combined_cluster_points = []

    for dbscan_label, coordinates in dbscan_cluster_dict.items():
        combined_cluster_points.extend(coordinates)
        combined_cluster_colors.extend([clusters_palette[( dbscan_label * 5 )%len(clusters_palette)]   if dbscan_label != -1 else [0, 0, 0, 0.1]] * len(coordinates) )

    ptcloud_all_clusters         = pv.PolyData(combined_cluster_points)
    ptcloud_all_clusters["rgba"] = combined_cluster_colors

    plotter.add_mesh(ptcloud_all_clusters, scalars="rgba", rgb=True, show_scalar_bar=False)

    # ? Visualize selected cluster
    plotter.subplot(0,1)

    lc = selected_cluster
    print("Max vals in selected cluster:", [[np.min(lc[:,0]), np.max(lc[:,0])], [np.min(lc[:,1]), np.max(lc[:,1])],[np.min(lc[:,2]), np.max(lc[:,2])] ])

    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster       = [[15, 10, 221, 1] for datapoint in selected_cluster]
    rgbas_positive      = np.array([[205, 209, 228, 0.2] for _ in positive_space])
    combined            = np.concatenate([selected_cluster, positive_space])
    rgbas_combined      = np.concatenate([rgbas_cluster, rgbas_positive])

    # list(selected_cluster).extend(list(positive_space))
    # list(rgbas_cluster).extend(list(rgbas_positive))
    # combined         = selected_cluster.extend(positive_space)
    # rgbas_combined   = np.concatenate([rgbas_cluster, rgbas_positie])
    container_points = [*list(selected_cluster),*list(positive_space)]
    # container_rgbas  = [*list(rgbas_cluster),*list(rgbas_positive)]

    np.save("selected_cluster.npy", selected_cluster)
    np.save("positive_space.npy", positive_space)
    point_cloud         = pv.PolyData(selected_cluster)
    point_cloud["rgba"] = rgbas_cluster

    point_cloud         = pv.PolyData(container_points)
    # point_cloud["rgba"] = container_rgbas
    plotter.add_points(point_cloud, rgb=True, show_scalar_bar=True)
    if gif:
        output_gif = gif_name
        # plotter.camera.zoom(1.5)
        plotter.open_gif(output_gif)

        # Rotate the camera 360 degrees
        for angle in range(0, 360, 2):  # 5 degree steps
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        plotter.show()

def visualize_mesh(mesh_path, rcsb_id: str | None = None, 
                  gif: bool = False, gif_name: str | None = None,
                  output_path: str | None = None):
    """
    Visualize a mesh file and optionally save to file.
    
    Parameters:
        mesh_path: Path to the mesh file
        rcsb_id: Optional structure ID for labeling
        gif: Whether to create a GIF
        gif_name: Name of the GIF file
        output_path: Path to save the visualization image
    """
    plotter = pv.Plotter(off_screen=output_path is not None or gif)

    _ = plotter.add_mesh(pv.read(mesh_path), opacity=0.8)
    plotter.add_axes(line_width=2, cone_radius=0.7, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.2, 0.8))
    plotter.add_text('RCSB_ID:{}'.format(rcsb_id if rcsb_id is not None else ""), 
                    position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    plotter.show_grid(n_xlabels=8, n_ylabels=8, n_zlabels=8, font_size=8)

    if output_path:
        plotter.screenshot(output_path)
        print(f"Saved visualization to {output_path}")

    if gif:
        output_gif = gif_name
        plotter.open_gif(output_gif)
        for angle in range(0, 360, 2):
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        if not output_path:  # Only show interactive window if not saving to file
            plotter.show()
    
    return plotter

def visualize_pointcloud(ptcloud, rcsb_id: str | None = None, 
                        gif: bool = False, gif_name: str | None = None,
                        output_path: str | None = None):
    """
    Visualize a point cloud and optionally save to file.
    
    Parameters:
        ptcloud: Point cloud data
        rcsb_id: Optional structure ID for labeling
        gif: Whether to create a GIF
        gif_name: Name of the GIF file
        output_path: Path to save the visualization image
    """
    plotter = pv.Plotter(off_screen=output_path is not None or gif)
    n_labels = 7
    plotter.show_grid(n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size=8)
    plotter.add_axes(line_width=2, cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.add_text('RCSB_ID:{}'.format(rcsb_id if rcsb_id is not None else ""), 
                    position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    plotter.add_points(pv.PolyData(ptcloud), color='black', point_size=3, opacity=0.3)
    
    if output_path:
        plotter.screenshot(output_path)
        print(f"Saved visualization to {output_path}")
        
    if gif:
        output_gif = gif_name
        plotter.open_gif(output_gif)
        for angle in range(0, 360, 2):
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        if not output_path:  # Only show interactive window if not saving to file
            plotter.show()
    
    return plotter

def visualize_pointcloud_axis(ptcloud,aux_ptcloud, base_point, axis_point, radius=0.1, height=None, rcsb_id:str|None=None, gif:bool=False, gif_name:str|None=None):
    # Convert inputs to numpy arrays to ensure compatibility
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)
    ptcloud = np.asarray(ptcloud)
    
    # Calculate the original axis vector
    original_axis = axis_point - base_point
    original_axis_length = np.linalg.norm(original_axis)
    original_axis = original_axis / original_axis_length
    
    # Set height to axis length if not specified
    if height is None:
        height = original_axis_length
    
    # Define the target z-axis (pointing down)
    target_z_axis = np.array([0, 0, -1])
    
    # Calculate the rotation matrix to align the original axis with the z-axis
    rotation_axis = np.cross(original_axis, target_z_axis)
    
    # Calculate the angle between the vectors
    angle = np.arccos(np.dot(original_axis, target_z_axis))
    
    # Create rotation matrix using Rodrigues' rotation formula
    if np.linalg.norm(rotation_axis) > 1e-10:
        rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
        K = np.array([
            [0, -rotation_axis[2], rotation_axis[1]],
            [rotation_axis[2], 0, -rotation_axis[0]],
            [-rotation_axis[1], rotation_axis[0], 0]
        ])
        R = (np.eye(3) + np.sin(angle) * K + 
             (1 - np.cos(angle)) * np.dot(K, K))
    else:
        # If the vectors are already aligned, use identity matrix
        R = np.eye(3)
    
    # Transform the point cloud
    transformed_ptcloud = np.dot(ptcloud - base_point, R.T)
    transformed_aux_pts = np.dot(aux_ptcloud - base_point, R.T)
    
    # Transform the base and axis points
    transformed_base_point = np.dot(base_point - base_point, R.T)
    transformed_axis_point = np.dot(axis_point - base_point, R.T)
    
    # Create the plotter
    plotter = pv.Plotter()
    
    # Add visualization elements
    n_labels = 7
    # plotter.show_grid(n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size=8)
    # plotter.add_axes(line_width=2, cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    
    # # Add RCSB ID text
    # plotter.add_text(f'RCSB_ID:{rcsb_id if rcsb_id is not None else ""}', 
    #                  position='upper_right', font_size=14, 
    #                  shadow=True, font='courier', color='black')
    
    # Create and add transformed cylinder
    cylinder = pv.Cylinder(
        center=transformed_base_point + (height/2) * np.array([0, 0, -1]),
        direction=[0, 0, -1],
        radius=radius,
        height=height,
        resolution=50
    ).triangulate()
    
    # Add cylinder with visible edges
    # plotter.add_mesh(cylinder, color='blue', opacity=0.3, show_edges=True, edge_color='black')

    plotter.add_mesh( cylinder, style='surface', line_width=3, show_edges=True, opacity=0.05, color='lightgreen', silhouette=True )
    plotter.add_points(transformed_base_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter.add_points(transformed_axis_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter.add_points(pv.PolyData(transformed_ptcloud), render_points_as_spheres=True,  style='points', point_size=3, opacity=0.4)
    plotter.add_points(pv.PolyData(transformed_aux_pts), render_points_as_spheres=True,  style='points', color='blue', point_size=5, opacity=0.8)
    
    # Show the plotr
    plotter.show()
    
    camera_position = plotter.camera_position
    # camera_focal_point = plotter.camera_focal_point
    # camera_viewup = plotter.camera_viewup

    # Then in a separate render
    plotter_render = pv.Plotter(off_screen=True)
    plotter_render.add_mesh( cylinder, style='surface', line_width=3, show_edges=True, opacity=0.05, color='lightgreen', silhouette=True )
    plotter_render.add_points(transformed_base_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter_render.add_points(transformed_axis_point, color='red', point_size=12, render_points_as_spheres=True)
    plotter_render.add_points(pv.PolyData(transformed_ptcloud), render_points_as_spheres=True,  style='points', point_size=3, opacity=0.4)
    plotter_render.add_points(pv.PolyData(transformed_aux_pts), render_points_as_spheres=True,  style='points', color='blue', point_size=5, opacity=0.8)
    
    plotter_render.camera_position = camera_position
    plotter_render.screenshot('4ug0_ptcloud.png', transparent_background=True, window_size=(1920,1280))
    
    return transformed_ptcloud

def create_cylinder_from_points(base_point, axis_point, radius, height):
    """
    Create a cylinder using two points: a base point and a point defining the axis direction.
    
    Parameters:
    base_point: array-like, [x, y, z] starting point of cylinder
    axis_point: array-like, [x, y, z] point defining cylinder axis direction
    radius: float, radius of cylinder
    height: float, total height of cylinder (can extend past axis_point)
    
    Returns:
    cylinder: PyVista cylinder mesh (triangulated)
    direction: numpy array of normalized direction vector
    """
    base_point = np.array(base_point)
    axis_point = np.array(axis_point)
    
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    ).triangulate()
    
    return cylinder, direction

def visualize_pointclouds(ptcloud1:np.ndarray, ptcloud2:np.ndarray, background_positive:np.ndarray, gif:bool=False, gif_name:str|None=None):
    plotter               = pv.Plotter(shape=(1, 2), off_screen=gif)
    plotter.subplot(0,0)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster1       = [[15, 100, 21, 1] for datapoint in ptcloud1]
    rgbas_positive1      = np.array([[205, 209, 228, 0.2] for _ in background_positive])

    combined1            = np.concatenate([ptcloud1, background_positive])
    rgbas_combined1      = np.concatenate([rgbas_cluster1, rgbas_positive1])

    point_cloud1         = pv.PolyData(combined1)
    point_cloud1["rgba"] = rgbas_combined1

    plotter.add_points(point_cloud1, scalars="rgba", rgb=True, show_scalar_bar=False)


    # ? Visualize selected cluster
    plotter.subplot(0,1)
    n_labels = 7
    plotter.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=1, ambient=1, label_size=(0.2, 0.6))
    plotter.show_grid( n_xlabels=n_labels, n_ylabels=n_labels, n_zlabels=n_labels, font_size = 12)

    rgbas_cluster2       = [[15, 10, 221, 1] for datapoint in ptcloud2]
    rgbas_positive2      = np.array([[205, 160, 200, 0.2] for _ in background_positive])
    combined2            = np.concatenate([ptcloud2, background_positive])
    rgbas_combined2      = np.concatenate([rgbas_cluster2, rgbas_positive2])
    point_cloud2         = pv.PolyData(combined2)
    point_cloud2["rgba"] = rgbas_combined2
    plotter.add_points(point_cloud2, scalars="rgba", rgb=True, show_scalar_bar=False)


    if gif:
        output_gif = gif_name
        # plotter.camera.zoom(1.5)
        plotter.open_gif(output_gif)

        # Rotate the camera 360 degrees
        for angle in range(0, 360, 2):  # 5 degree steps
            plotter.camera.azimuth = angle
            plotter.write_frame()
        plotter.close()
        print(f"GIF saved as {output_gif}")
    else:
        plotter.show()


    plotter.show()

def visualize_clipping_result(
    original_points: np.ndarray,
    clipped_points: np.ndarray,
    mesh_path: str,
    show_mesh: bool = True,
):
    """
    Visualizes the original points, clipped points, and the clipping mesh.

    Parameters:
        original_points (np.ndarray): The original point cloud
        clipped_points (np.ndarray): The clipped point cloud
        mesh_path (str): Path to the mesh file
        show_mesh (bool): Whether to show the mesh or not
    """
    p = pv.Plotter()

    # Add original points in red
    original_cloud = pv.PolyData(original_points)
    p.add_mesh(
        original_cloud,
        color="red",
        point_size=5,
        render_points_as_spheres=True,
        label="Original Points",
    )

    # Add clipped points in blue
    clipped_cloud = pv.PolyData(clipped_points)
    p.add_mesh(
        clipped_cloud,
        color="blue",
        point_size=5,
        render_points_as_spheres=True,
        label="Clipped Points",
    )

    # Add mesh if requested
    if show_mesh:
        mesh = pv.read(mesh_path)
        p.add_mesh(
            mesh, style="wireframe", color="gray", opacity=0.5, label="Clipping Mesh"
        )

    p.add_legend()
    p.show()

def visualize_filtered_residues(
    filtered_residues: List[T],
    all_residues: Optional[List[T]],
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray] = lambda x: x.center_of_mass(),
    point_size: float = 5,
    opacity: float = 0.3,
    show: bool = True,
    output_path: Optional[str] = None,
    window_size: tuple = (1024, 768)
) -> Optional[pv.Plotter]:
    """
    Visualize filtered residues alongside the cylinder that was used for filtering.
    """
    # Initialize plotter
    plotter = pv.Plotter(window_size=window_size, off_screen=output_path is not None and not show)
    
    # Get positions of filtered residues
    filtered_positions = np.array([
        position_getter(res) for res in filtered_residues
    ])
    
    # Create and add cylinder
    direction = axis_point - base_point
    direction = direction / np.linalg.norm(direction)
    
    cylinder = pv.Cylinder(
        center=base_point + (height/2) * direction,
        direction=direction,
        radius=radius,
        height=height,
        resolution=30
    )
    
    # Add cylinder with wireframe and solid style
    plotter.add_mesh(
        cylinder, 
        style='wireframe', 
        color='red', 
        line_width=2, 
        label='Cylinder'
    )
    plotter.add_mesh(
        cylinder, 
        style='surface', 
        color='red', 
        opacity=opacity,
    )
    
    # If all_residues provided, show unfiltered residues in gray
    if all_residues is not None:
        all_positions = np.array([
            position_getter(res) for res in all_residues
        ])
        plotter.add_points(
            all_positions, 
            color='gray', 
            point_size=point_size-2, 
            opacity=0.3,
            label='Unfiltered Residues'
        )
    
    # Add filtered points
    plotter.add_points(
        filtered_positions, 
        color='blue', 
        point_size=point_size,
        label='Filtered Residues'
    )
    
    # Add base and axis points for reference
    plotter.add_points(
        np.array([base_point]), 
        color='green', 
        point_size=point_size*2, 
        label='Base Point'
    )
    plotter.add_points(
        np.array([axis_point]), 
        color='yellow', 
        point_size=point_size*2, 
        label='Axis Point'
    )
    
    # Add axis line
    axis_line = pv.Line(base_point, axis_point)
    plotter.add_mesh(
        axis_line, 
        color='white', 
        line_width=2, 
        label='Cylinder Axis'
    )
    
    # Add legend
    plotter.add_legend()
    
    # Set camera position for better initial view
    plotter.camera_position = 'iso'
    plotter.enable_eye_dome_lighting()  # Improves point visibility
    
    # Add text with stats
    stats_text = (
        f'Total filtered residues: {len(filtered_residues)}\n'
        f'Cylinder height: {height:.1f}\n'
        f'Cylinder radius: {radius:.1f}'
    )
    if all_residues:
        stats_text = f'Total residues: {len(all_residues)}\n' + stats_text
        
    plotter.add_text(
        stats_text,
        position='upper_left',
        font_size=12,
        color='white'
    )
    
    # Save screenshot if output path is provided
    if output_path:
        plotter.screenshot(output_path)
        print(f"Saved visualization to {output_path}")
    
    if show:
        plotter.show()
        return None
        
    return plotter

def ptcloud_convex_hull_points_and_gif(
    pointcloud: np.ndarray,
    ALPHA: float,
    TOLERANCE: float,
    output_path: str,
    n_frames: int = 180,
) -> np.ndarray:

    from tqdm import tqdm
    from PIL import Image

    assert pointcloud is not None

    # Create the convex hull
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(alpha=ALPHA, tol=TOLERANCE, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface()

    # Set up the plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(convex_hull, color="lightblue", show_edges=True)
    plotter.add_points(cloud, color="red", point_size=5)

    # Set up the camera
    plotter.camera_position = "xy"

    # Create frames
    frames = []
    for i in tqdm(range(n_frames)):
        plotter.camera.azimuth = i * (360 / n_frames)
        plotter.render()
        image = plotter.screenshot(transparent_background=False, return_img=True)
        frames.append(Image.fromarray(image))

    # Save as GIF
    frames[0].save(
        output_path, save_all=True, append_images=frames[1:], duration=50, loop=0
    )

    print(f"GIF saved to {output_path}")

    return convex_hull.points

def estimate_normals_and_create_gif(
    convex_hull_surface_pts: np.ndarray,
    output_path: str,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
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

    frames[0].save(
        output_path, save_all=True, append_images=frames[1:], duration=50, loop=0
    )

    print(f"GIF saved to {output_path}")
