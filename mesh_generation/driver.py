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
import numpy as np
from sklearn.cluster import DBSCAN
from __archive.scripts.pymol_visualtion import extract_chains
from mesh_generation.bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from compas.geometry import bounding_box
from mesh_generation.visualization import DBSCAN_CLUSTERS_visualize_largest, custom_cluster_recon_path, plot_multiple_by_kingdom, plot_multiple_surfaces, plot_with_landmarks, DBSCAN_CLUSTERS_particular_eps_minnbrs
from mesh_generation.paths import *
from mesh_generation.voxelize import (expand_atomcenters_to_spheres_threadpool, normalize_atom_coordinates)
from ribctl import EXIT_TUNNEL_WORK, POISSON_RECON_BIN, RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.libpdb import extract_chains_by_auth_asym_id


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

def index_grid(expanded_sphere_voxels: np.ndarray):

    normalized_sphere_cords, mean_abs_vectors = normalize_atom_coordinates(expanded_sphere_voxels)
    voxel_size = 1
    sphere_cords_quantized = np.round( np.array(normalized_sphere_cords / voxel_size) ).astype(int)
    max_values      = np.max(sphere_cords_quantized, axis=0)
    grid_dimensions = max_values + 1
    vox_grid        = np.zeros(grid_dimensions)

    vox_grid[
        sphere_cords_quantized[:, 0],
        sphere_cords_quantized[:, 1],
        sphere_cords_quantized[:, 2],
    ] = 1

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

def surface_pts_via_convex_hull(
    rcsb_id: str, 
    selected_cluster: np.ndarray 
)->np.ndarray:
    assert selected_cluster is not None
    cloud       = pv.PolyData(selected_cluster)
    grid        = cloud.delaunay_3d(alpha=2, tol=1.5, offset=2, progress_bar=True)
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points


def save_mesh_point_cloud_as_ply( mesh, save_path: str):
    pcd        = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(mesh)
    o3d.io.write_point_cloud(save_path, pcd)
    print("Wrote {}".format(save_path))


def estimate_normals(rcsb_id, convex_hull_surface_pts: np.ndarray):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamHybrid(radius=5, max_nn=15) )
    pcd.orient_normals_consistent_tangent_plane(k=10)
    o3d.io.write_point_cloud(surface_with_normals_path(rcsb_id), pcd)

    print("Wrote surface with normals {}".format(surface_with_normals_path(rcsb_id)))

def pick_largest_poisson_cluster(clusters_container:dict[int,list])->np.ndarray:
    DBSCAN_CLUSTER_ID = 1
    for k, v in clusters_container.items():
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)

        # print("Picked cluster {} because it has more points({})".format(DBSCAN_CLUSTER_ID, len(clusters_container[DBSCAN_CLUSTER_ID])))
    return np.array(clusters_container[DBSCAN_CLUSTER_ID])

def apply_poisson_reconstruction(rcsb_id: str, poisson_recon_path:str):
    import plyfile

    command   = [ POISSON_RECON_BIN, "--in", surface_with_normals_path(rcsb_id), "--out", poisson_recon_path, "--depth", "6", "--pointWeight", "3", ]
    process   = subprocess.run(command, capture_output=True, text=True)

    if process.returncode == 0:
        print("PoissonRecon executed successfully.")
        print("WRote {}".format(poisson_recon_path))
        # Convert the plyfile to asciii
        data      = plyfile.PlyData.read(poisson_recon_path)
        data.text = True
        data.write(poisson_recon_path + ".txt")
        print("Wrote {}".format(poisson_recon_path + ".txt"))
    else:
        print("Error:", process.stderr)


#TODO
def save_lsu_alpha_chains(rcsb_id:str, reconstructed_tunnel_ply:str, outpath:str)->str:
    
    # grab all the chains that are within the NeighborSearch of the tunnel 5 angstrom
    # extract them from the structure with pymol, save to disc

    import pyvista as pv
    import numpy as np

    # Load the PLY file using PyVista
    mesh = pv.read(reconstructed_tunnel_ply)

    # Extract points as a NumPy array
    points_array = np.array(mesh.points)

    mmcif_parser = MMCIFParser(QUIET=True)
    structure    = mmcif_parser.get_structure(rcsb_id, os.path.join(RIBETL_DATA,rcsb_id,rcsb_id + ".cif"))
    atoms        = list(structure.get_atoms())
    ns           = NeighborSearch(atoms)

    # Define a distance threshold for the neighbor search
    distance_threshold = 5.0

    neighbor_chains_auth_asym_ids = set()
    # Perform the neighbor search
    for point in points_array:
        _ = ns.search(point, distance_threshold)
        [neighbor_chains_auth_asym_ids.add(chain_name) for chain_name in [ a.get_full_id()[2] for a in _]]



    extract_chains_by_auth_asym_id(rcsb_id,list( neighbor_chains_auth_asym_ids ), outpath)
    return outpath


def construct_trimming_alphashape(rcsb_id:str, lsu_chains_file:str, alpha,tol):
    mmcif_parser = MMCIFParser(QUIET=True)
    structure    = mmcif_parser.get_structure(rcsb_id+"alpha", lsu_chains_file)
    atoms        = np.array([a.get_coord() for a in list(structure.get_atoms())])
    cloud        = pv.PolyData(atoms)
    delaunay_shape         = cloud.delaunay_3d(alpha=alpha, tol=tol, offset=2, progress_bar=True)
    surface_polydata = delaunay_shape.extract_surface()
    convex_hull = delaunay_shape.extract_surface().cast_to_pointset()
    return delaunay_shape,convex_hull.points, surface_polydata

#TODO
def trim_with_alphasurface(rcsb_id:str,reconstruction_pcl:np.ndarray, alpha:float)->np.ndarray:
    """
    Trim the poisson reconstruction with the alpha shape surface
    """
    pass


def ____pipeline(RCSB_ID):
    _u_EPSILON     = 5.5
    _u_MIN_SAMPLES = 600
    _u_METRIC      = "euclidean"

    struct_tunnel_dir = os.path.join(EXIT_TUNNEL_WORK, RCSB_ID)
    if not os.path.exists(struct_tunnel_dir):
        os.mkdir(struct_tunnel_dir)

    if not os.path.exists(spheres_expanded_pointset_path(RCSB_ID)):
        "the data arrives here as atom coordinates extracted from the biopython model "
        if not os.path.exists(tunnel_atom_encoding_path(RCSB_ID)):
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
        print("Saved expanded sphere pointset to {}".format(spheres_expanded_pointset_path(RCSB_ID)))

    else:
        bbox_atoms_expanded = np.load(spheres_expanded_pointset_path(RCSB_ID))

    xyz_positive, xyz_negative, _ , translation_vectors = index_grid(bbox_atoms_expanded) 
    np.save(translation_vectors_path(RCSB_ID), translation_vectors)
    db, clusters_container = interior_capture_DBSCAN( xyz_negative, _u_EPSILON, _u_MIN_SAMPLES, _u_METRIC )
    largest_cluster = pick_largest_poisson_cluster(clusters_container)

    #! Transform the cluster back into original coordinate frame
    coordinates_in_the_original_frame =  largest_cluster  - translation_vectors[1] + translation_vectors[0]
    main_cluster = pv.PolyData(coordinates_in_the_original_frame)

    
    delaunay, ptcloud, surface_polydata = construct_trimming_alphashape(RCSB_ID, mmcif_ensemble_LSU(RCSB_ID), alpha=8, tol=3)
    selected = main_cluster.select_enclosed_points(surface_polydata)
    pts      = main_cluster.extract_points(selected['SelectedPoints'].view(bool),adjacent_cells=False)

    print(pts)
    print(np.shape(pts.points))

    pl = pv.Plotter()
    _ = pl.add_mesh(main_cluster, style='wireframe')
    _ = pl.add_mesh(LSU_alpha_shape, style='wireframe')
    _ = pl.add_points(pts, color='r')
    pl.show()

    exit(1)

    surface_pts     = surface_pts_via_convex_hull( RCSB_ID, coordinates_in_the_original_frame )
    np.save(convex_hull_cluster_path(RCSB_ID), surface_pts)
    estimate_normals(RCSB_ID, surface_pts)
    apply_poisson_reconstruction(RCSB_ID, poisson_recon_path(RCSB_ID))
        
def main():

    # ? ---------- Params ------------
    parser = argparse.ArgumentParser()
    # Add command-line arguments
    parser.add_argument( "--rcsb_id", type=str, help="Specify the value for eps (float)", required=True )
    parser.add_argument( "--eps", type=float, help="Specify the value for eps (float)")
    parser.add_argument( "--min_samples", type=int, help="Specify the value for min_samples (int)" )
    parser.add_argument( "--metric", choices=DBSCAN_METRICS, help="Choose a metric from the provided options", )

    parser.add_argument( "--full_pipeline",   action='store_true')
    parser.add_argument( "--final",   action='store_true')
    parser.add_argument( "--dbscan",   action='store_true')
    parser.add_argument( "--dbscan_tuple",  type=str)

    parser.add_argument( "--multisurf",   action='store_true')
    parser.add_argument( "--fig",   action='store_true')
    parser.add_argument( "--kingdom",   choices=['bacteria','archaea','eukaryota'])

    parser.add_argument( "--trim",   action='store_true')

    args          = parser.parse_args()
    RCSB_ID       = args.rcsb_id.upper()


    if args.trim:
        import pyvista as pv

        outpath           = '{}/{}/{}_lsu_alphashape.mmcif'.format(EXIT_TUNNEL_WORK, RCSB_ID, RCSB_ID)
        alpha_chains_path = '{}/{}/{}_lsu_alphashape.ply'.format(EXIT_TUNNEL_WORK, RCSB_ID, RCSB_ID)

        if not os.path.exists(outpath):
            save_lsu_alpha_chains(RCSB_ID, poisson_recon_path(RCSB_ID), outpath)
        else:
            ...
           
        delaunay, point_cloud = construct_trimming_alphashape(RCSB_ID, outpath, alpha=8, tol=3)
        delaunay.save(alpha_shape_LSU(RCSB_ID))
        poisson_recon_mesh    = pv.read(poisson_recon_path(RCSB_ID))

        # plotter     = pv.Plotter(shape=(1,2))
        # plotter.subplot(0,0)
        # plotter.add_mesh(poisson_recon_mesh, opacity=1)
        # plotter.add_mesh(delaunay, opacity=0.4, color="blue")
        # # plotter.add_mesh(point_cloud, opacity=0.5, color="red")
        # plotter.subplot(0,1)
        # plotter.add_mesh(delaunay, opacity=1, color="blue")
        # plotter.show()

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
        ____pipeline(RCSB_ID)

    if args.final:

        eps,min_nbrs =  args.dbscan_tuple.split(",")
        eps = float(eps)
        min_nbrs = int(min_nbrs)
        plot_with_landmarks(RCSB_ID, float(eps),int(min_nbrs),custom_cluster_recon_path(RCSB_ID, float(eps), int(min_nbrs)))

    if args.multisurf:
        plot_multiple_surfaces(RCSB_ID)

    if args.kingdom:
        eps,min_nbrs =  args.dbscan_tuple.split(",")
        eps = float(eps)
        min_nbrs = int(min_nbrs)
        plot_multiple_by_kingdom(args.kingdom, eps, min_nbrs)
if __name__ == "__main__":
    main()
