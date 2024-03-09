
# ? ---------- Paths ------------
import os
from ribctl import EXIT_TUNNEL_WORK, RIBETL_DATA


tunnel_atom_encoding_path = lambda rcsb_id: os.path.join(
    EXIT_TUNNEL_WORK,
    rcsb_id.upper(),
    "{}_tunnel_atoms_bbox.json".format(rcsb_id.upper()),
)
spheres_expanded_pointset_path = lambda rcsb_id: os.path.join(
    EXIT_TUNNEL_WORK,
    rcsb_id.upper(),
    "{}_spheres_expanded_pointset.npy".format(rcsb_id.upper()),
)
translation_vectors_path = lambda rcsb_id: os.path.join(
    EXIT_TUNNEL_WORK,
    rcsb_id.upper(),
    "{}_normalization_vectors.npy".format(rcsb_id.upper()),
)
selected_dbscan_cluster_path = lambda rcsb_id: os.path.join( EXIT_TUNNEL_WORK, rcsb_id.upper(), "{}_dbscan_cluster.npy".format(rcsb_id.upper()) )
convex_hull_cluster_path     = lambda rcsb_id: os.path.join( EXIT_TUNNEL_WORK, rcsb_id.upper(), "{}_convex_hull.npy".format(rcsb_id.upper()) )
surface_with_normals_path    = lambda rcsb_id: os.path.join( EXIT_TUNNEL_WORK, rcsb_id.upper(), "{}_normal_estimated_surf.ply".format(rcsb_id.upper()), )
poisson_recon_path           = lambda rcsb_id: os.path.join( EXIT_TUNNEL_WORK, rcsb_id.upper(), "{}_poisson_recon.ply".format(rcsb_id.upper()) )
ptc_data_path                = lambda rcsb_id: os.path.join( RIBETL_DATA, rcsb_id.upper(), "{}_PTC_COORDINATES.json".format(rcsb_id.upper()) )

def mmcif_ensemble_LSU(rcsb_id):
    return os.path.join(EXIT_TUNNEL_WORK, rcsb_id, '{}_lsu_alphashape.mmcif'.format(rcsb_id))

def convex_hull_ensemble_LSU(rcsb_id):
    return os.path.join(EXIT_TUNNEL_WORK, rcsb_id, '{}_lsu_convex_hull.ply'.format(rcsb_id))

def alpha_shape_LSU(rcsb_id):
    return os.path.join(EXIT_TUNNEL_WORK, rcsb_id, '{}_lsu_alphashape.ply'.format(rcsb_id))

def custom_cluster_recon_path(rcsb_id, eps, min_nbrs):
    return os.path.join( EXIT_TUNNEL_WORK, rcsb_id.upper(), "{}_poisson_recon-eps{}_minnbrs{}.ply".format(rcsb_id.upper(), eps, min_nbrs) )
# ? ------------------------------