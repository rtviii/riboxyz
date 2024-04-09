from pprint import pprint
from matplotlib import pyplot as plt
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
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, DBSCAN_CLUSTERS_particular_eps_minnbrs, visualize_mesh, visualize_pointcloud
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, index_grid)
from ribctl import EXIT_TUNNEL_WORK
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_lsu_ensemble_tunnel_vicinity



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

    centerline_expansion_atoms       = parse_struct_via_centerline( rcsb_id, open_tunnel_csv(rcsb_id) )
    centerline_expansion_coordinates = np.array( [a.get_coord() for a in centerline_expansion_atoms] )

    bbox                             = bounding_box(centerline_expansion_coordinates)
    bbox_interior_atoms              = parse_struct_via_bbox(rcsb_id, bbox)

    return encode_atoms(
        rcsb_id,
        bbox_interior_atoms,
        write     = True,
        writepath = tunnel_atom_encoding_path(rcsb_id) )

def expand_bbox_atoms_to_spheres(atom_coordinates:np.ndarray, sphere_vdw_radii:np.ndarray, rcsb_id: str):
    sphere_sources = zip(atom_coordinates, sphere_vdw_radii)
    SINK           = []
    expanded       = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)

    return np.array(expanded)


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
    return np.array(clusters_container[DBSCAN_CLUSTER_ID])



def pipeline(RCSB_ID,args):

    #! [ Pipeline Parameters ]
    
    _u_EPSILON     = 5.5 if args.dbscan_tuple is None else float(args.dbscan_tuple.split(",")[0])
    _u_MIN_SAMPLES = 600 if args.dbscan_tuple is None else int(args.dbscan_tuple.split(",")[1])
    _u_METRIC      = "euclidean"

    d3d_alpha      = args.D3D_alpha if args.D3D_alpha is not None else 2
    d3d_tol        = args.D3D_tol if args.D3D_tol is not None else 1
    PR_depth       = args.PR_depth if args.PR_depth is not None else 6
    PR_ptweight    = args.PR_ptweight if args.PR_ptweight is not None else 3

    struct_tunnel_dir = os.path.join(EXIT_TUNNEL_WORK, RCSB_ID)

    if not os.path.exists(struct_tunnel_dir):
        os.mkdir(struct_tunnel_dir)

    #! [ Bounding Box Atoms ]
    if args.bbox or ( not os.path.exists(spheres_expanded_pointset_path(RCSB_ID)) ) :
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
    else:
        print("Opened atoms in the bounding box")
        bbox_atoms_expanded = np.load(spheres_expanded_pointset_path(RCSB_ID))




    # ! [ Bounding Box Atoms are transformed into an Index Grid ]
    # _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded)
    voxel_grid, grid_dimensions, translation_vectors = index_grid(bbox_atoms_expanded)


    #! [ Extract the largest cluster from the DBSCAN clustering ]
    db, clusters_container = interior_capture_DBSCAN(np.asarray(np.where(voxel_grid != 1)).T, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
    largest_cluster        = pick_largest_poisson_cluster(clusters_container)

    #! [ Visualize the largest DBSCAN cluster to establish whether trimming is required ]
    visualize_pointcloud(largest_cluster, RCSB_ID)



    if args.trim:
        # user_input = input("Truncate bbox? Enter tuples of the format 'x,20 : z,40 : y,15 : Y,20' (lowercase for truncation from origin, uppercase for truncation from end of axis) or 'Q' to quit: ")
        user_input = input("Truncate bbox? Enter tuples of the format ' 10:69|20:80|5:70' (for x|y|z axis truncation, ||20:50 to skip axis) or 'Q' to quit: ")
        if user_input.lower() == 'q':
            print("Exiting the program.")
            exit(0)
        try:
            truncation_strings = user_input.replace(" ", '').split("|")
            if len(truncation_strings) != 3:
                print("You have to enter three truncation parameters for x, y, and z axes. (ex. ||20:50 or to skip x and y axes,  or |40:-| to cut y from 40 to the 'end' np style)")

            x_tuple = [ int(c) if c != '' else None for c  in truncation_strings[0].split(":") ] if ":" in truncation_strings[0] else None
            y_tuple = [ int(c) if c != '' else None for c  in truncation_strings[1].split(":") ] if ":" in truncation_strings[1] else None
            z_tuple = [ int(c) if c != '' else None for c  in truncation_strings[2].split(":") ] if ":" in truncation_strings[2] else None
            
        except Exception as e:
            print("Failed to parse trim parameters:" , e)
            exit(1)

    TRUNCATION_TUPLES = [x_tuple, y_tuple, z_tuple]
    print("Got truncation_tuples:", TRUNCATION_TUPLES)

    if TRUNCATION_TUPLES is not None:
        if len(TRUNCATION_TUPLES) != 3:
            raise IndexError("You have to enter three truncation parameters for x, y, and z axes. (ex. ||20:50 to skip x and y axes )")

        if TRUNCATION_TUPLES[0] is not None:
            if TRUNCATION_TUPLES[0][0] is not None:
                voxel_grid[:TRUNCATION_TUPLES[0][0],:,:]  = 1

            if TRUNCATION_TUPLES[0][1] is not None:
                 voxel_grid[TRUNCATION_TUPLES[0][1]:,:,:] = 1


        if TRUNCATION_TUPLES[1] is not None:
            if TRUNCATION_TUPLES[1][0] is not None:
                voxel_grid[:,:TRUNCATION_TUPLES[1][0],:]  = 1

            if TRUNCATION_TUPLES[1][1] is not None:
                 voxel_grid[:,TRUNCATION_TUPLES[1][1]:,:] = 1



        if TRUNCATION_TUPLES[2] is not None:

            if TRUNCATION_TUPLES[2][0] is not None:
                voxel_grid[:,:,:TRUNCATION_TUPLES[2][0]]  = 1

            if TRUNCATION_TUPLES[2][1] is not None:
                 voxel_grid[:,:,TRUNCATION_TUPLES[2][1]:] = 1

    # xyz_positive, xyz_negative, grid_dimensions , translation_vectors = index_grid(bbox_atoms_expanded, TRUNCATION_TUPLES=[x_tuple, y_tuple, z_tuple] if args.trim else None)
    print("Truncated pointcloud")
    visualize_pointcloud(np.asarray(np.where(voxel_grid != 1)).T, RCSB_ID)
    visualize_pointcloud(np.asarray(np.where(voxel_grid == 1)).T, RCSB_ID)

    db, clusters_container = interior_capture_DBSCAN( np.asarray(np.where(voxel_grid != 1)).T, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
    largest_cluster = pick_largest_poisson_cluster(clusters_container)
    DBSCAN_CLUSTERS_visualize_largest(np.asarray(np.where(voxel_grid == 1)).T, clusters_container, largest_cluster)

    #! [ Transform the cluster back into Original Coordinate Frame ]
    coordinates_in_the_original_frame = largest_cluster  - translation_vectors[1] + translation_vectors[0]
    surface_pts                       = ptcloud_convex_hull_points(coordinates_in_the_original_frame, 3,2)
    np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
    apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=6, recon_pt_weight=3)
    visualize_mesh(pv.read(poisson_recon_path(RCSB_ID)), RCSB_ID)
