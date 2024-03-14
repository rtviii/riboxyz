import argparse
from pprint import pprint
import subprocess
import typing
from matplotlib import pyplot as plt
import open3d as o3d
from Bio.PDB import MMCIFParser
from Bio.PDB import NeighborSearch
import pyvista as pv
import json
import os
import numpy as np
from sklearn.cluster import DBSCAN
from __archive.scripts.pymol_visualtion import extract_chains
from mesh_generation.bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from compas.geometry import bounding_box
from mesh_generation.libsurf import apply_poisson_reconstruction, estimate_normals, ptcloud_convex_hull_points
from mesh_generation.lsu_alpha_surface import lsu_ensemble_convex_hull, lsu_ensemble_get_chains, ptcloud_convex_hull, vestibule_sphere_expansion
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, DBSCAN_CLUSTERS_particular_eps_minnbrs
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, normalize_atom_coordinates)
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_lsu_ensemble_tunnel_vicinity


DBSCAN_METRICS        = [
    "braycurtis",
    "canberra",
    "chebyshev",
    "correlation",
    "dice",
    "hamming",
    "jaccard",
    "kulsinski",
    "mahalanobis", # ?
    "minkowski",   # ?
    "rogerstanimoto",
    "russellrao",
    "seuclidean",
    "sokalmichener",
    "sokalsneath",
    "sqeuclidean", #?
    "yule",
    "cityblock",
    "cosine",
    "euclidean",
    "l1",
    "l2",
    "manhattan",
]

def bounding_box(points: np.ndarray):
    """Computes the axis-aligned minimum bounding box of a list of points.

    Parameters
    ----------
    points : sequence[point]
        XYZ coordinates of the points.

    Returns
    -------
    list[[float, float, float]]
        XYZ coordinates of 8 points defining a box.


    """
    x, y, z = zip(*points)
    min_x = min(x)
    max_x = max(x)
    min_y = min(y)
    max_y = max(y)
    min_z = min(z)
    max_z = max(z)
    return [
        [min_x, min_y, min_z],
        [max_x, min_y, min_z],
        [max_x, max_y, min_z],
        [min_x, max_y, min_z],
        [min_x, min_y, max_z],
        [max_x, min_y, max_z],
        [max_x, max_y, max_z],
        [min_x, max_y, max_z],
    ]

def extract_bbox_atoms(rcsb_id: str) -> list:

    centerline_expansion_atoms = parse_struct_via_centerline(
        rcsb_id, open_tunnel_csv(rcsb_id)
    )
    centerline_expansion_coordinates = np.array(
        [a.get_coord() for a in centerline_expansion_atoms]
    )
    bbox                = bounding_box(centerline_expansion_coordinates)
    bbox_interior_atoms = parse_struct_via_bbox(rcsb_id, bbox)

    return encode_atoms(
        rcsb_id,
        bbox_interior_atoms,
        write=True,
        writepath=tunnel_atom_encoding_path(rcsb_id),
    )

def expand_bbox_atoms_to_spheres(atom_coordinates:np.ndarray, sphere_vdw_radii:np.ndarray, rcsb_id: str):

    sphere_sources = zip(atom_coordinates, sphere_vdw_radii)
    SINK = []
    expanded = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)

    return np.array(expanded)

def index_grid(expanded_sphere_voxels: np.ndarray, TRUNCATION_TUPLES:list[tuple[str, int]]|None = None) :

    normalized_sphere_cords, mean_abs_vectors = normalize_atom_coordinates(expanded_sphere_voxels)
    voxel_size = 1

    sphere_cords_quantized = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)
    max_values             = np.max(sphere_cords_quantized, axis=0)
    grid_dimensions        = max_values + 1
    vox_grid               = np.zeros(grid_dimensions)

    vox_grid[
        sphere_cords_quantized[:, 0],
        sphere_cords_quantized[:, 1],
        sphere_cords_quantized[:, 2],
    ] = 1


    if TRUNCATION_TUPLES is not None:
        for ( TRUNCATION_AXIS, TRUNCATION_FACTOR) in TRUNCATION_TUPLES:
            match TRUNCATION_AXIS: 
                case "x":
                    vox_grid = vox_grid[:-TRUNCATION_FACTOR,:,:]
                case "X":
                    vox_grid = vox_grid[TRUNCATION_FACTOR:,:,:]
                case "y":
                    vox_grid = vox_grid[:,:-TRUNCATION_FACTOR,:]
                case "Y":
                    vox_grid = vox_grid[:,TRUNCATION_FACTOR:,:]
                case "z":
                    vox_grid = vox_grid[:,:,:-TRUNCATION_FACTOR]
                case "Z":
                    vox_grid = vox_grid[:,:,TRUNCATION_FACTOR:]
                case _:
                    print("Invalid truncation axis. Please use 'x', 'X', 'y', 'Y', 'z', 'Z'.")
        print("Truncated {} by {}".format(TRUNCATION_AXIS, TRUNCATION_FACTOR))

    print("Voxel grid shape(post truncation if applied): ", vox_grid.shape)
    __xyz_v_negative_ix = np.asarray(np.where(vox_grid != 1))
    __xyz_v_positive_ix = np.asarray(np.where(vox_grid == 1))


    return (
        __xyz_v_positive_ix.T,
        __xyz_v_negative_ix.T,
        grid_dimensions,
        mean_abs_vectors,
    )

def interior_capture_DBSCAN(
    xyz_v_negative: np.ndarray,
    eps           ,
    min_samples   ,
    metric        : str = "euclidean",
): 

    cluster_colors = dict(zip(range(-1, 40), plt.cm.terrain(np.linspace(0, 1, 40))))
    for k, v in cluster_colors.items():
        cluster_colors[k] = [*v[:3], 0.5]

    u_EPSILON     = eps
    u_MIN_SAMPLES = min_samples
    u_METRIC      = metric

    print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(xyz_v_negative), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) ) 

    db     = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit( xyz_v_negative )
    labels = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(xyz_v_negative, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))

    return db, CLUSTERS_CONTAINER

def pick_largest_poisson_cluster(clusters_container:dict[int,list])->np.ndarray:
    DBSCAN_CLUSTER_ID = 1
    for k, v in clusters_container.items():
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)

        # print("Picked cluster {} because it has more points({})".format(DBSCAN_CLUSTER_ID, len(clusters_container[DBSCAN_CLUSTER_ID])))
    return np.array(clusters_container[DBSCAN_CLUSTER_ID])

def ____pipeline(RCSB_ID,args):

    _u_EPSILON     = 5.5 if args.dbscan_tuple is None else float(args.dbscan_tuple.split(",")[0])
    _u_MIN_SAMPLES = 600 if args.dbscan_tuple is None else int(args.dbscan_tuple.split(",")[1])
    _u_METRIC      = "euclidean"

    struct_tunnel_dir = os.path.join(EXIT_TUNNEL_WORK, RCSB_ID)
    if not os.path.exists(struct_tunnel_dir):
        os.mkdir(struct_tunnel_dir)

    if not os.path.exists(spheres_expanded_pointset_path(RCSB_ID)):
        "the data arrives here as atom coordinates extracted from the biopython model "
        if not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)) and args.bbox_radius != None:
            bbox_atoms          = extract_bbox_atoms(RCSB_ID)

            _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
            _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))
            bbox_atoms_expanded = expand_bbox_atoms_to_spheres(_atom_centers,_vdw_radii, RCSB_ID)

        else:
            with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile: 
                bbox_atoms = json.load(infile)
            _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
            _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))
            bbox_atoms_expanded = expand_bbox_atoms_to_spheres(_atom_centers, _vdw_radii, RCSB_ID)

        np.save(spheres_expanded_pointset_path(RCSB_ID), bbox_atoms_expanded)
        print(">>Saved expanded sphere pointset to {}".format(spheres_expanded_pointset_path(RCSB_ID)))

    else:
        bbox_atoms_expanded = np.load(spheres_expanded_pointset_path(RCSB_ID))


    # ! index grid

    _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded)
    # np.save(translation_vectors_path(RCSB_ID), translation_vectors)

    #! truncate the bounding box:
    db, clusters_container = interior_capture_DBSCAN( xyz_negative, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
    largest_cluster        = pick_largest_poisson_cluster(clusters_container)

    pl =  pv.Plotter()
    pl.add_points(largest_cluster, color='r', point_size=2, render_points_as_spheres=True)
    pl.add_axes(line_width=2,cone_radius=0.3, shaft_length=2, tip_length=0.3, ambient=0.5, label_size=(0.4, 0.16),)
    pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    pl.show_grid()
    pl.show()



    #! Visualize the cluster to establish whether trimming is required


    if args.trim:
        user_input = input("Truncate bbox? Enter tuples of the format 'x,20:z,40:Y,20' or 'Q' to quit: ")
        if user_input.lower() == 'q':
            print("Exiting the program.")
            exit(0)
        try:

            truncation_string = user_input.split(":")
            truncation_params = [( str( pair.split(",")[0] ) , int( pair.split(",")[1] )) for pair in truncation_string]

            _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded,TRUNCATION_TUPLES=truncation_params) 
            # np.save(translation_vectors_path(RCSB_ID), translation_vectors)

            #! truncate the bounding box:
            # pl =  pv.Plotter()
            # pl.add_points(xyz_negative, color='r', point_size=2)
            # pl.show()

            db, clusters_container = interior_capture_DBSCAN( xyz_negative, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
            largest_cluster = pick_largest_poisson_cluster(clusters_container)

            #! Transform the cluster back into original coordinate frame
            coordinates_in_the_original_frame =  largest_cluster  - translation_vectors[1] + translation_vectors[0]

            # #! Vestibule truncation via surface
            # main_cluster              = pv.PolyData(coordinates_in_the_original_frame)
            # vestibule_expansion_mesh_ = pv.read(alphashape_ensemble_LSU(RCSB_ID))
            # selected                  = main_cluster.select_enclosed_points(vestibule_expansion_mesh_, check_surface=True)
            # pts                       = main_cluster.extract_points( selected['SelectedPoints'].view(bool), adjacent_cells=False, )
            # exit()

            surface_pts     = ptcloud_convex_hull_points(coordinates_in_the_original_frame, 2.5,1)
            print(">>Extracted convex hull points.")
            print(np.array(surface_pts).shape)
            np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
            estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
            print(">>Estimated normals")
            apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=6, recon_pt_weight=3)
            
            pl                        = pv.Plotter()
            _                         = pl.add_mesh(pv.read(poisson_recon_path(RCSB_ID)), opacity=0.8)
            pl.add_axes(line_width=5,cone_radius=0.6, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.4, 0.16),)
            pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
            pl.show_grid()
            pl.show()
            print(f"You truncation parameters: {stepsize}, {axis}")

        except Exception as e:
            print("Invalid input. Please enter an integer or 'Q' to quit.: ", e)


    #! Transform the cluster back into original coordinate frame
    coordinates_in_the_original_frame =  largest_cluster  - translation_vectors[1] + translation_vectors[0]


    # #! Vestibule truncation via surface
    # main_cluster              = pv.PolyData(coordinates_in_the_original_frame)
    # vestibule_expansion_mesh_ = pv.read(alphashape_ensemble_LSU(RCSB_ID))
    # selected                  = main_cluster.select_enclosed_points(vestibule_expansion_mesh_, check_surface=True)
    # pts                       = main_cluster.extract_points( selected['SelectedPoints'].view(bool), adjacent_cells=False, )
    # exit()

    surface_pts     = ptcloud_convex_hull_points(coordinates_in_the_original_frame, 3,2)
    print(">>Extracted convex hull points.")
    print(np.array(surface_pts).shape)
    np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
    print(">>Estimated normals")
    apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=6, recon_pt_weight=3)
    
    pl                        = pv.Plotter()
    _                         = pl.add_mesh(pv.read(poisson_recon_path(RCSB_ID)), opacity=0.8)
    pl.add_axes(line_width=5,cone_radius=0.6, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.4, 0.16),)
    pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
    pl.show_grid()
    pl.show()

        

        
def main():

    # ? ---------- Params ------------
    parser = argparse.ArgumentParser()
    # Add command-line arguments

    parser.add_argument( "--full_pipeline",   action='store_true')

    # visualization options
    parser.add_argument( "--final",   action='store_true')
    parser.add_argument( "--dbscan",   action='store_true')
    parser.add_argument( "--multisurf",   action='store_true')
    parser.add_argument( "--kingdom",   choices=['bacteria','archaea','eukaryota'])

        # truncation
    parser.add_argument( "--lsu_alpha",   action='store_true')


    # pipeline parameters
    parser.add_argument( "--rcsb_id", type=str, help="Specify the value for eps (float)", required=True )
    parser.add_argument( "--metric", choices=DBSCAN_METRICS, help="Choose a metric from the provided options", )
    parser.add_argument( "--bbox_radius",  type=int, help="The radius of the bbox expansion", required=False)
    parser.add_argument( "--trim",  action='store_true',required=False)
    parser.add_argument( "--truncation_factor", type=int, help="How many voxels (A) to leave trim on the vestibule side.")
    parser.add_argument( "--dbscan_tuple",  type=str)


    args          = parser.parse_args()
    RCSB_ID       = args.rcsb_id.upper()


    if args.lsu_alpha:


        ALPHA_VAL = 7.9
        ALPHA_TOL = 2
       
        vestibule_sphere_ptcloud = vestibule_sphere_expansion(RCSB_ID, 50)
        convex_hull = ptcloud_convex_hull(vestibule_sphere_ptcloud, ALPHA_VAL, ALPHA_TOL, offset=2)

        pcd = estimate_normals(convex_hull.points, convex_hull_ensemble_LSU(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15,correction_tangent_planes_n=15)

        # pcd        = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(convex_hull.points))
        o3d.visualization.draw_geometries([pcd])
        radii=[2    ,8]
        rec_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting( pcd, o3d.utility.DoubleVector(radii))
        o3d.visualization.draw_geometries([pcd, rec_mesh])


        apply_poisson_reconstruction(convex_hull_ensemble_LSU(RCSB_ID), alphashape_ensemble_LSU(RCSB_ID), recon_depth=6, recon_pt_weight=7)


        FONT                  = 'courier'
        RCSB_ID="6Z6K"
        plotter               = pv.Plotter()
        mesh_  = pv.read(alphashape_ensemble_LSU(RCSB_ID))
        plotter.add_mesh(mesh_, opacity=0.5)
        plotter.add_text('ALPHA VAL: {}'.format(ALPHA_VAL), position='upper_left', font_size=20, shadow=True, font=FONT, color='black')
        plotter.show()

    if args.dbscan:
        if args.dbscan_tuple is not None:
            eps,min_nbrs       =  args.dbscan_tuple.split(",")
            metric= 'euclidean'
            # expand_bbox_atoms_to_spheres(atom_coordinates:np.ndarray, sphere_vdw_radii:np.ndarray, rcsb_id: str):
            expanded_sphere_voxels = np.load(spheres_expanded_pointset_path(RCSB_ID))
            xyz_pos, xyz_neg, _, _= index_grid(expanded_sphere_voxels)
            db,clusters_container = interior_capture_DBSCAN( xyz_neg,  float(eps), int(min_nbrs), metric)

            largest_cluster = pick_largest_poisson_cluster(clusters_container)
            surface_pts     = surface_pts_via_convex_hull( RCSB_ID, largest_cluster )
            np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
            apply_poisson_reconstruction(RCSB_ID, custom_cluster_recon_path(RCSB_ID, eps, min_nbrs))
            # DBSCAN_CLUSTERS_visualize_largest(xyz_pos, clusters_container, largest_cluster)
            # DBSCAN_CLUSTERS_visualize_largest(xyz_pos, clusters_container, largest_cluster)
            DBSCAN_CLUSTERS_particular_eps_minnbrs(clusters_container, float(eps),int(min_nbrs))

    if args.full_pipeline:
        ____pipeline(RCSB_ID, args)

    if args.final:

        eps,min_nbrs =  args.dbscan_tuple.split(",")
        eps = float(eps)
        min_nbrs = int(min_nbrs)
        plot_with_landmarks(RCSB_ID, float(eps),int(min_nbrs),poisson_recon_path(RCSB_ID))

    if args.multisurf:
        plot_multiple_surfaces(RCSB_ID)

    if args.kingdom:
        eps,min_nbrs =  args.dbscan_tuple.split(",")
        eps = float(eps)
        min_nbrs = int(min_nbrs)
        plot_multiple_by_kingdom(args.kingdom, eps, min_nbrs)
if __name__ == "__main__":
    main()
