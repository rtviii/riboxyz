import os
import numpy as np
import pyvista as pv
import open3d as o3d

# from __scripts.cylinder import transform_points_to_C0
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.kdtree_approach import (
    create_point_cloud_mask,
    transform_points_from_C0,
    transform_points_to_C0,
)
from ribctl.lib.npet.mesh_full_pipeline import (
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)
from ribctl.lib.npet.mesh_libsurf import (
    apply_poisson_reconstruction,
    estimate_normals,
    ptcloud_convex_hull_points,
)
from ribctl.lib.npet.mesh_visualization import (
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
    visualize_mesh,
    visualize_pointcloud,
)
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager
from ribctl.lib.npet.tunnel_bbox_ptc_constriction import (
    filter_residues_parallel,
    ribosome_entities,
)
from ribctl.lib.npet.util import landmark_constriction_site, landmark_ptc


def create_npet_mesh(RCSB_ID: str):
    assets                 = TunnelMeshAssetsManager(RCSB_ID)
    cifpath                = AssetType.MMCIF.get_path(RCSB_ID)

    ashapepath             = AssetType.ALPHA_SHAPE.get_path(RCSB_ID)
    meshpath               = AssetType.NPET_MESH.get_path(RCSB_ID)

    print(f"Debug - meshpath: {meshpath}, type: {type(meshpath)}")
    R                      = 30
    H                      = 100
    Vsize                  = 1
    ATOM_SIZE              = 2
    normals_pcd_path       = assets.tunnel_pcd_normal_estimated
    ashape_watertight_mesh = pv.read(ashapepath)

    if not os.path.exists(ashapepath):
        raise FileNotFoundError(f"File {ashapepath} not found")
    # # All atoms
    ptc_pt             = np.array(landmark_ptc(RCSB_ID))
    constriction_pt    = np.array(landmark_constriction_site(RCSB_ID))
    entities           = ribosome_entities(RCSB_ID, cifpath, "R")
    filtered           = filter_residues_parallel(entities, ptc_pt, constriction_pt, R, H)
    filtered_points    = np.array( [atom.get_coord() for residue in filtered for atom in residue.child_list] )
    transformed_points = transform_points_to_C0( filtered_points, ptc_pt, constriction_pt )

    mask, (x, y, z) = create_point_cloud_mask( transformed_points, radius=R, height=H, voxel_size=Vsize, radius_around_point=ATOM_SIZE, )

    points            = np.where(~mask)
    empty_coordinates = np.column_stack((x[points[0]], y[points[1]], z[points[2]]))

    back_projected        = transform_points_from_C0( empty_coordinates, ptc_pt, constriction_pt )
    select                = pv.PolyData(back_projected).select_enclosed_points(ashape_watertight_mesh)
    mask                  = select["SelectedPoints"]
    interior              = back_projected[mask == 1]
    empty_in_world_coords = np.array(interior)

    _u_EPSILON_initial_pass     = 5.5
    _u_MIN_SAMPLES_initial_pass = 600
    _u_EPSILON_refinement       = 3.5
    _u_MIN_SAMPLES_refinement   = 175

    d3d_alpha   = 2
    d3d_tol     = 1
    PR_depth    = 6
    PR_ptweight = 3

    db, clusters_container = DBSCAN_capture( empty_in_world_coords, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass ) #! [ Extract the largest cluster from the DBSCAN clustering ]
    largest_cluster, largest_cluster_id = DBSCAN_pick_largest_cluster(
        clusters_container
    )
    # visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs( clusters_container, _u_EPSILON_initial_pass, _u_MIN_SAMPLES_initial_pass, ptc_pt, constriction_pt, np.array([[1,1,1]]))
    # exit()
    db, refined_clusters_container = DBSCAN_capture(
        largest_cluster, _u_EPSILON_refinement, _u_MIN_SAMPLES_refinement
    )
    refined_cluster, refined_cluster_id = DBSCAN_pick_largest_cluster(
        refined_clusters_container
    )

    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
        clusters_container,
        _u_EPSILON_initial_pass,
        _u_MIN_SAMPLES_initial_pass,
        ptc_pt,
        constriction_pt,
        refined_cluster,
        R,
        H,
    )

    #! [ Transform the cluster back into original coordinate frame ]
    surface_pts = ptcloud_convex_hull_points(refined_cluster, d3d_alpha, d3d_tol)
    visualize_pointcloud(surface_pts, RCSB_ID)

    #! [ Transform the cluster back into Original Coordinate Frame ]
    normal_estimated_pcd = estimate_normals(
        surface_pts,
        kdtree_radius=10,
        kdtree_max_nn=15,
        correction_tangent_planes_n=10,
    )

    o3d.io.write_point_cloud(normals_pcd_path, normal_estimated_pcd)
    print("Wrote to {}".format(normals_pcd_path))

    apply_poisson_reconstruction(
        normals_pcd_path,
        meshpath,
        recon_depth=PR_depth,
        recon_pt_weight=PR_ptweight,
    )

    visualize_mesh(meshpath)
