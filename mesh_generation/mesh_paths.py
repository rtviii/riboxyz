# ? ---------- Paths ------------
import os
from ribctl import EXIT_TUNNEL_WORK, RIBETL_DATA
from ribctl.etl.assets_structure import StructureAssets

# most of this stuff is unnecessary
tunnel_atom_encoding_path      = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_tunnel_atoms_bbox.json".format(rcsb_id.upper()) )
spheres_expanded_pointset_path = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_spheres_expanded_pointset.npy".format(rcsb_id.upper()), )
# translation_vectors_path       = lambda rcsb_id: os.path.join( Assets(rcsb_id).paths.tunnel_dir, "{}_normalization_vectors.npy".format(rcsb_id.upper()), )
selected_dbscan_cluster_path   = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_dbscan_cluster.npy".format(rcsb_id.upper()) )
convex_hull_cluster_path       = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_convex_hull.npy".format(rcsb_id.upper()) )
surface_with_normals_path      = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_normal_estimated_surf.ply".format(rcsb_id.upper()), )
poisson_recon_path             = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_poisson_recon.ply".format(rcsb_id.upper()) )
poisson_recon_ascii_path       = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.tunnel_dir, "{}_poisson_recon_ascii.ply".format(rcsb_id.upper()) )
ptc_data_path                  = lambda rcsb_id: os.path.join( StructureAssets(rcsb_id).paths.ptc)

TRIMMING_PARAMS_DICT_PATH = os.path.join( EXIT_TUNNEL_WORK, "trimming_params_dict.json" )
TUNNEL_PATH               = lambda rcsb_id :  os.path.join( EXIT_TUNNEL_WORK, "mole_tunnels", "tunnel_{}.csv".format(rcsb_id) )

def custom_cluster_recon_path(rcsb_id, eps, min_nbrs):
    return os.path.join( StructureAssets(rcsb_id).paths.ptc, rcsb_id.upper(), "{}_poisson_recon-eps{}_minnbrs{}.ply".format(rcsb_id.upper(), eps, min_nbrs) )
