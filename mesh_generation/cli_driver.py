import argparse
from pprint import pprint
from matplotlib import pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
from mesh_generation.bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from mesh_generation.paths import *
from mesh_generation.full_pipeline import pipeline


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

# def ____pipeline(RCSB_ID,args):

#     _u_EPSILON     = 5.5 if args.dbscan_tuple is None else float(args.dbscan_tuple.split(",")[0])
#     _u_MIN_SAMPLES = 600 if args.dbscan_tuple is None else int(args.dbscan_tuple.split(",")[1])
#     _u_METRIC      = "euclidean"


#     d3d_alpha   = args.D3D_alpha if args.D3D_alpha is not None else 2
#     d3d_tol     = args.D3D_tol if args.D3D_tol is not None else 1
#     PR_depth    = args.PR_depth if args.PR_depth is not None else 6
#     PR_ptweight = args.PR_ptweight if args.PR_ptweight is not None else 3


#     struct_tunnel_dir = os.path.join(EXIT_TUNNEL_WORK, RCSB_ID)
#     if not os.path.exists(struct_tunnel_dir):
#         os.mkdir(struct_tunnel_dir)

#     if not os.path.exists(spheres_expanded_pointset_path(RCSB_ID)):
#         "the data arrives here as atom coordinates extracted from the biopython model "
#         if not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)) and args.bbox_radius != None:
#             bbox_atoms          = extract_bbox_atoms(RCSB_ID)

#             _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
#             _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))
#             bbox_atoms_expanded = expand_bbox_atoms_to_spheres(_atom_centers,_vdw_radii, RCSB_ID)

#         else:
#             with open( tunnel_atom_encoding_path(RCSB_ID), "r", ) as infile: 
#                 bbox_atoms = json.load(infile)
#             _atom_centers       = np.array(list(map(lambda x: x["coord"], bbox_atoms)))
#             _vdw_radii          = np.array(list(map(lambda x: x["vdw_radius"], bbox_atoms)))
#             bbox_atoms_expanded = expand_bbox_atoms_to_spheres(_atom_centers, _vdw_radii, RCSB_ID)

#         np.save(spheres_expanded_pointset_path(RCSB_ID), bbox_atoms_expanded)
#         print(">>Saved expanded sphere pointset to {}".format(spheres_expanded_pointset_path(RCSB_ID)))

#     else:
#         bbox_atoms_expanded = np.load(spheres_expanded_pointset_path(RCSB_ID))


#     # ! index grid

#     _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded)
#     # np.save(translation_vectors_path(RCSB_ID), translation_vectors)

#     #! truncate the bounding box:
#     db, clusters_container = interior_capture_DBSCAN( xyz_negative, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
#     largest_cluster        = pick_largest_poisson_cluster(clusters_container)

#     pl =  pv.Plotter()
#     pl.add_points(largest_cluster, color='r', point_size=2, render_points_as_spheres=True)
#     pl.add_axes(line_width=4,cone_radius=0.7, shaft_length=2, tip_length=0.9, ambient=0.5, label_size=(0.4, 0.16),)
#     pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
#     pl.show_grid(
#          n_xlabels=6,
#         n_ylabels=6,
#         n_zlabels=6,
#     )
   
#     pl.show()



#     #! Visualize the cluster to establish whether trimming is required


#     if args.trim:
#         user_input = input("Truncate bbox? Enter tuples of the format 'x,20:z,40:Y,20' or 'Q' to quit: ")
#         if user_input.lower() == 'q':
#             print("Exiting the program.")
#             exit(0)
#         try:

#             truncation_string = user_input.split(":")
#             truncation_params = [( str( pair.split(",")[0] ) , int( pair.split(",")[1] )) for pair in truncation_string]

#             _, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded,TRUNCATION_TUPLES=truncation_params) 
#             np.save(translation_vectors_path(RCSB_ID), translation_vectors)

#             db, clusters_container = interior_capture_DBSCAN( xyz_negative, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
#             largest_cluster = pick_largest_poisson_cluster(clusters_container)

#             #! Transform the cluster back into original coordinate frame
#             coordinates_in_the_original_frame =  largest_cluster  - translation_vectors[1] + translation_vectors[0]

#             surface_pts     = ptcloud_convex_hull_points(coordinates_in_the_original_frame, d3d_alpha ,d3d_tol)
#             print(">>Extracted convex hull points.")
#             print(np.array(surface_pts).shape)
#             np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
#             estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
#             print(">>Estimated normals")
#             apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=PR_depth, recon_pt_weight=PR_ptweight)
            
#             pl                        = pv.Plotter()
#             _                         = pl.add_mesh(pv.read(poisson_recon_path(RCSB_ID)), opacity=0.8)
#             pl.add_axes(line_width=5,cone_radius=0.6, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.4, 0.16),)
#             pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
#             pl.show_grid(
#                 n_xlabels = 10,
#                 n_ylabels = 10,
#                 n_zlabels = 10,
#             )
#             pl.show()
#             exit(0)

#         except Exception as e:
#             print("Invalid input. Please enter an integer or 'Q' to quit.: ", e)


#     #! Transform the cluster back into original coordinate frame
#     coordinates_in_the_original_frame =  largest_cluster  - translation_vectors[1] + translation_vectors[0]


#     # #! Vestibule truncation via surface
#     # main_cluster              = pv.PolyData(coordinates_in_the_original_frame)
#     # vestibule_expansion_mesh_ = pv.read(alphashape_ensemble_LSU(RCSB_ID))
#     # selected                  = main_cluster.select_enclosed_points(vestibule_expansion_mesh_, check_surface=True)
#     # pts                       = main_cluster.extract_points( selected['SelectedPoints'].view(bool), adjacent_cells=False, )
#     # exit()

#     surface_pts     = ptcloud_convex_hull_points(coordinates_in_the_original_frame, 3,2)
#     print(">>Extracted convex hull points.")
#     print(np.array(surface_pts).shape)
#     np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
#     estimate_normals(surface_pts, surface_with_normals_path(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15, correction_tangent_planes_n=10)
#     print(">>Estimated normals")
#     apply_poisson_reconstruction(surface_with_normals_path(RCSB_ID), poisson_recon_path(RCSB_ID), recon_depth=6, recon_pt_weight=3)
    
#     pl                        = pv.Plotter()
#     _                         = pl.add_mesh(pv.read(poisson_recon_path(RCSB_ID)), opacity=0.8)
#     pl.add_axes(line_width=5,cone_radius=0.6, shaft_length=0.7, tip_length=0.3, ambient=0.5, label_size=(0.4, 0.16),)
#     pl.add_text('RCSB_ID:{}'.format(RCSB_ID), position='upper_right', font_size=14, shadow=True, font='courier', color='black')
#     pl.show_bounds( n_xlabels=6, n_ylabels=6, n_zlabels=6, )
#     pl.show()

        

        
def main():

    # ? ---------- Params ------------
    parser = argparse.ArgumentParser()
    # Add command-line arguments

    parser.add_argument( "--full_pipeline",   action='store_true')

    # visualization options
    # parser.add_argument( "--final",   action='store_true')
    # parser.add_argument( "--dbscan",   action='store_true')
    # parser.add_argument( "--multisurf",   action='store_true')
    # parser.add_argument( "--kingdom",   choices=['bacteria','archaea','eukaryota'])

    parser.add_argument( "--metric", choices=DBSCAN_METRICS, help="Choose a metric from the provided options", )
    # truncation
    # parser.add_argument( "--lsu_alpha",   action='store_true')

    # pipeline parameters
    parser.add_argument( "--rcsb_id", type=str, help="Specify the value for eps (float)", required=True )
    parser.add_argument( "--bbox",  action='store_true', help="Extract the bounding box atoms and save them to a file")
    parser.add_argument( "--bbox_radius",  type=int, help="The radius of the bbox expansion", required=False)

    #! Reconstruction Parameters
    #?      - dbsan
    parser.add_argument( "--dbscan_tuple",  type=str)

    #?      - delaunay_3d
    parser.add_argument( "--D3D_alpha",  type=float)
    parser.add_argument( "--D3D_tol",  type=float)

    #?      - poisson Reconstruction
    parser.add_argument( "--PR_depth",  type=int)
    parser.add_argument( "--PR_ptweight",  type=int)

    #! Interactive / Visualization
    parser.add_argument( "--trim",  action='store_true',required=False)

    args          = parser.parse_args()
    RCSB_ID       = args.rcsb_id.upper()

    if args.full_pipeline:
        pipeline(RCSB_ID, args)
        exit(0)


if __name__ == "__main__":
    main()


    # if args.dbscan:
    #     if args.dbscan_tuple is not None:
    #         eps,min_nbrs       =  args.dbscan_tuple.split(",")
    #         metric= 'euclidean'
    #         # expand_bbox_atoms_to_spheres(atom_coordinates:np.ndarray, sphere_vdw_radii:np.ndarray, rcsb_id: str):
    #         expanded_sphere_voxels = np.load(spheres_expanded_pointset_path(RCSB_ID))
    #         xyz_pos, xyz_neg, _, _= index_grid(expanded_sphere_voxels)
    #         db,clusters_container = interior_capture_DBSCAN( xyz_neg,  float(eps), int(min_nbrs), metric)

    #         largest_cluster = pick_largest_poisson_cluster(clusters_container)
    #         surface_pts     = surface_pts_via_convex_hull( RCSB_ID, largest_cluster )
    #         np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    #         apply_poisson_reconstruction(RCSB_ID, custom_cluster_recon_path(RCSB_ID, eps, min_nbrs))
    #         # DBSCAN_CLUSTERS_visualize_largest(xyz_pos, clusters_container, largest_cluster)
    #         # DBSCAN_CLUSTERS_visualize_largest(xyz_pos, clusters_container, largest_cluster)
    #         DBSCAN_CLUSTERS_particular_eps_minnbrs(clusters_container, float(eps),int(min_nbrs))

    # if args.multisurf:
    #     plot_multiple_surfaces(RCSB_ID)

    # if args.kingdom:
    #     eps,min_nbrs =  args.dbscan_tuple.split(",")
    #     eps = float(eps)
    #     min_nbrs = int(min_nbrs)
    #     plot_multiple_by_kingdom(args.kingdom, eps, min_nbrs)




    # if args.lsu_alpha:

    #     ALPHA_VAL = 7.9
    #     ALPHA_TOL = 2
       
    #     vestibule_sphere_ptcloud = vestibule_sphere_expansion(RCSB_ID, 50)
    #     convex_hull              = ptcloud_convex_hull(vestibule_sphere_ptcloud, ALPHA_VAL, ALPHA_TOL, offset=2)
    #     pcd                      = estimate_normals(convex_hull.points, convex_hull_ensemble_LSU(RCSB_ID), kdtree_radius=10, kdtree_max_nn=15,correction_tangent_planes_n=15)

    #     # pcd        = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(convex_hull.points))
    #     # o3d.visualization.draw_geometries([pcd])
    #     # radii=[2    ,8]
    #     # rec_mesh = o3d.geometry.TriangleMesh.create_from_point_cloud_ball_pivoting( pcd, o3d.utility.DoubleVector(radii))
    #     # o3d.visualization.draw_geometries([pcd, rec_mesh])


    #     apply_poisson_reconstruction(convex_hull_ensemble_LSU(RCSB_ID), alphashape_ensemble_LSU(RCSB_ID), recon_depth=6, recon_pt_weight=7)

    #     FONT                  = 'courier'
    #     RCSB_ID="6Z6K"
    #     plotter               = pv.Plotter()
    #     mesh_  = pv.read(alphashape_ensemble_LSU(RCSB_ID))
    #     plotter.add_mesh(mesh_, opacity=0.5)
    #     plotter.add_text('ALPHA VAL: {}'.format(ALPHA_VAL), position='upper_left', font_size=20, shadow=True, font=FONT, color='black')
    #     plotter.show()
