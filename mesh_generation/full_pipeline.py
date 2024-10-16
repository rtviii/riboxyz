import pyvista as pv
import json
import os
import numpy as np
from sklearn.cluster import DBSCAN
from mesh_generation.bbox_extraction import ( encode_atoms, extract_bbox_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from mesh_generation.libsurf import apply_poisson_reconstruction, estimate_normals, estimate_normals_and_create_gif, ptcloud_convex_hull_points, ptcloud_convex_hull_points_and_gif
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs, visualize_mesh, visualize_pointcloud, visualize_pointclouds 
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, index_grid)
import numpy as np

def expand_bbox_atoms_to_spheres(atom_coordinates:np.ndarray, sphere_vdw_radii:np.ndarray, rcsb_id: str):
    sphere_sources = zip(atom_coordinates, sphere_vdw_radii)
    SINK           = []
    expanded       = expand_atomcenters_to_spheres_threadpool(SINK, sphere_sources)

    return np.array(expanded)

def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps           ,
    min_samples   ,
    metric        : str = "euclidean",
): 

    u_EPSILON     = eps
    u_MIN_SAMPLES = min_samples
    u_METRIC      = metric

    print( "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format( len(ptcloud), u_EPSILON, u_MIN_SAMPLES, u_METRIC ) ) 

    db     = DBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit( ptcloud )
    labels = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(ptcloud, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))
    return db, CLUSTERS_CONTAINER

def DBSCAN_pick_largest_cluster(clusters_container:dict[int,list], pick_manually:bool=False)->np.ndarray:
    DBSCAN_CLUSTER_ID = 0
    if pick_manually:
        print("Running Manual Cluster")
        picked_id =  int(input("Enter Cluster ID to proceed the reconstruction with\n (options:[{}]):".format(list( clusters_container.keys() ))))
        return np.array(clusters_container[picked_id])

    for k, v in clusters_container.items():
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)
    return np.array(clusters_container[DBSCAN_CLUSTER_ID])

def cache_trimming_parameters( RCSB_ID:str, trim_tuple:list, file_path=TRIMMING_PARAMS_DICT_PATH):
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File {file_path} not found.")

    with open(file_path, 'r') as file:
        data = json.load(file)

    data[RCSB_ID] = trim_tuple
    with open(file_path, 'w') as file:
        json.dump(data, file, indent=4)
        print(f"Entry '{RCSB_ID}' added to {file_path}")

def load_trimming_parameters( RCSB_ID:str, file_path=TRIMMING_PARAMS_DICT_PATH):
    with open(file_path, 'r') as file:
        data = json.load(file)
    return data[RCSB_ID] if RCSB_ID in data else None

def pipeline(RCSB_ID,args):
    GIF_INTERMEDIATES = args.gif_intermediates


    #! [ Pipeline Parameters ]
    _u_EPSILON     = 5.5 if args.dbscan_tuple is None else float(args.dbscan_tuple.split(",")[0])
    _u_MIN_SAMPLES = 600 if args.dbscan_tuple is None else int(args.dbscan_tuple.split(",")[1])
    _u_METRIC      = "euclidean"

    d3d_alpha         = args.D3D_alpha   if args.D3D_alpha   is not None else 2
    d3d_tol           = args.D3D_tol     if args.D3D_tol     is not None else 1
    PR_depth          = args.PR_depth    if args.PR_depth    is not None else 6
    PR_ptweight       = args.PR_ptweight if args.PR_ptweight is not None else 3

    struct_tunnel_dir = Assets(RCSB_ID).paths.tunnel_dir
    if not os.path.exists(struct_tunnel_dir):
        os.mkdir(struct_tunnel_dir)

    if args.bbox or  ( not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)) ):
        extract_bbox_atoms(RCSB_ID)

    #! [ Bounding Box Atoms ]
    if args.bbox or ( not os.path.exists(spheres_expanded_pointset_path(RCSB_ID)) ) :
        "the data arrives here as *bound_box around the centerline expansion* atom coordinates extracted from the biopython model "

        with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile: 
            bbox_atoms = json.load(infile)
        _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
        _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))
        bbox_atoms_expanded = expand_bbox_atoms_to_spheres(_atom_centers, _vdw_radii, RCSB_ID)

        np.save(spheres_expanded_pointset_path(RCSB_ID), bbox_atoms_expanded)
        print("Saved spheres_expanded_pointset data to disk : {}".format(spheres_expanded_pointset_path(RCSB_ID)))
    else:
        print("Opened atoms in the bounding box")
        bbox_atoms_expanded = np.load(spheres_expanded_pointset_path(RCSB_ID))

    
    visualize_pointcloud(bbox_atoms_expanded, RCSB_ID, GIF_INTERMEDIATES, "{}.bbox_atoms_expanded.gif".format(RCSB_ID))
    # ! [ Bounding Box Atoms are transformed into an Index Grid ]
    # _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded)
    initial_grid, grid_dimensions, translation_vectors = index_grid(bbox_atoms_expanded)


    # visualize_pointcloud(initial_grid, RCSB_ID, False, "{}.initial.grid.gif".format(RCSB_ID))
    # ? Here no trimming has yet occurred.

    #! [ Invert the grid ]
    inverted_grid = np.asarray(np.where(initial_grid != 1)).T
    
    # visualize_pointcloud(inverted_grid, RCSB_ID, False, "{}.initial.grid.gif".format(RCSB_ID))
    #! [ Capture DBSCAN clusters within the "negative" space]
    db, clusters_container = DBSCAN_capture(inverted_grid, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC ) 
    #! [ Extract the largest cluster from the DBSCAN clustering ]
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(clusters_container, _u_EPSILON, _u_MIN_SAMPLES, GIF_INTERMEDIATES, "{}.dbscan.clusters.gif".format(RCSB_ID))
    largest_cluster = DBSCAN_pick_largest_cluster(clusters_container)
    

    # #! [ Visualize the largest DBSCAN cluster to establish whether trimming is required ]
    # DBSCAN_CLUSTERS_visualize_largest(np.asarray(np.where(initial_grid == 1)).T, clusters_container, largest_cluster)
    visualize_pointcloud(largest_cluster, RCSB_ID, GIF_INTERMEDIATES, "{}.ptcloud.gif".format(RCSB_ID))

    # TODO : refactor this trimming logic out
    TRUNCATION_TUPLES = load_trimming_parameters(RCSB_ID)
    if args.trim and TRUNCATION_TUPLES is None:
        # user_input = input("Truncate bbox? Enter tuples of the format 'x,20 : z,40 : y,15 : Y,20' (lowercase for truncation from origin, uppercase for truncation from end of axis) or 'Q' to quit: ")
        user_input = input("Truncate bbox? Enter tuples of the format ' 10:69|20:80|5:70' (for x|y|z axis truncation, ||20:50 to skip axis) or 'Q' to quit: ")
        if user_input.lower() == 'q':
            print("Exiting the program.")
            exit(0)
        try:
            truncation_strings = user_input.replace(" ", '').split("|")

            if len(truncation_strings) != 3:
                print("You have to enter three truncation parameters for x, y, and z axes. (ex. ||20:50 or to skip x and y axes,  or |40:-| to cut y from 40 to the 'end' np style)")

            x_tuple = [ int(xx) if xx != '' else None for xx  in truncation_strings[0].split(":") ] if ":" in truncation_strings[0] else None
            y_tuple = [ int(yy) if yy != '' else None for yy  in truncation_strings[1].split(":") ] if ":" in truncation_strings[1] else None
            z_tuple = [ int(zz) if zz != '' else None for zz  in truncation_strings[2].split(":") ] if ":" in truncation_strings[2] else None
            
        except Exception as e:
            print("Failed to parse trim parameters:" , e)
            exit(1)

        TRUNCATION_TUPLES = [x_tuple, y_tuple, z_tuple]
        cache_trimming_parameters(RCSB_ID, TRUNCATION_TUPLES)
    if TRUNCATION_TUPLES is not None:
        if len(TRUNCATION_TUPLES) != 3:
            raise IndexError("You have to enter three truncation parameters for x, y, and z axes. (ex. ||20:50 to skip x and y axes )")

        def trim_pt_filter(pt:np.ndarray):

            if TRUNCATION_TUPLES[0] is not None:
                assert(len(TRUNCATION_TUPLES[0]) == 2)
                x_start,x_end = TRUNCATION_TUPLES[0]
                if x_start is not None:
                    if pt[0] <= x_start:
                        return False
                if x_end is not None:
                    if pt[0] >= x_end:
                        return False

            if TRUNCATION_TUPLES[1] is not None:
                assert(len(TRUNCATION_TUPLES[1]) == 2)
                y_start,y_end = TRUNCATION_TUPLES[1]

                if y_start is not None:
                    if pt[1] <= y_start:
                        return False

                if y_end is not None:
                    if pt[1] >= y_end:
                        return False

            if TRUNCATION_TUPLES[2] is not None:
                assert(len(TRUNCATION_TUPLES[2]) == 2)
                z_start,z_end = TRUNCATION_TUPLES[2]

                if z_start is not None:
                    if pt[2] <= z_start:
                        return False

                elif z_end is not None:
                    if pt[2] >= z_end:
                        return False

            return True
    if args.trim:
        trimmed_cluster = np.array(list(filter(trim_pt_filter,list(largest_cluster))))
    else:
        trimmed_cluster = largest_cluster
    # TODO : refactor this trimming logic out


    # Visualize
    # visualize_pointcloud(trimmed_cluster, RCSB_ID, GIF_INTERMEDIATES, "{}.ptcloud_trimmed.gif".format(RCSB_ID))
    # ! ----------
    # ! [ Extract the largest DBSCAN cluster with more restrictive parameters so as to capture the geometry more precisely and avoid trimmed the merging of trimmed parts into the main cluster   ]
    # ! Second pass of dbscan is needed to pick up the largest part of the trimmed cluster. 
    # ! (it's not uncommon that the truncated shape contains more than one disconnected cluster. Ex. imagine cutting a crescent moon in half in the coronal plane) 

    _, dbscan_container= DBSCAN_capture(trimmed_cluster, 3  , 123, _u_METRIC)

    print("DBSCAN Clusters: ")
    for (k,v) in dbscan_container.items():
        print(k, len(v))

    # visualize_pointcloud(main_cluster, RCSB_ID, GIF_INTERMEDIATES, "{}.ptcloud_trimmed_sharpened.gif".format(RCSB_ID))
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(dbscan_container, 3, 123, GIF_INTERMEDIATES, "{}.ptcloud_trimmed_sharpened.gif".format(RCSB_ID))
    main_cluster = DBSCAN_pick_largest_cluster(dbscan_container, args.cluster_manual)
    # ! ----------

    #! [ Transform the cluster back into Original Coordinate Frame ]
    coordinates_in_the_original_frame = main_cluster  - translation_vectors[1] + translation_vectors[0]

    #! [ Transform the cluster back into original coordinate frame ]
    surface_pts = ptcloud_convex_hull_points(coordinates_in_the_original_frame, d3d_alpha,d3d_tol)
    visualize_pointcloud(surface_pts, RCSB_ID, False, "{}.surface_pts.gif".format(RCSB_ID))

    #! [ Transform the cluster back into Original Coordinate Frame ]
    np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
    apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=PR_depth, recon_pt_weight=PR_ptweight)

    visualize_mesh(pv.read(poisson_recon_path(RCSB_ID)), RCSB_ID, GIF_INTERMEDIATES, "{}.reconstruction.gif".format(RCSB_ID))
