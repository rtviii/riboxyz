import warnings
from Bio.PDB.MMCIFParser import MMCIFParser
import numpy as np
import os
import alphashape
import trimesh

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager

data_dir = os.environ.get("DATA_DIR")
warnings.filterwarnings("ignore")


def cif_to_point_cloud(cif_path: str):
    parser    = MMCIFParser()
    structure = parser.get_structure("structure", cif_path)
    coordinates = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    coordinates.append(atom.coord)  # Get Cartesian coordinates
    return coordinates


def sample_within_alpha_shape(alpha_shape_mesh, num_samples=1000):
    """
    Samples points within the boundaries of a 3D alpha shape.

    Parameters:
    - alpha_shape_mesh: trimesh.Trimesh
        The alpha shape mesh that defines the 3D boundaries.
    - num_samples: int
        The number of points to sample within the alpha shape.

    Returns:
    - sampled_points: np.ndarray
        Array of sampled points within the alpha shape of shape (num_samples, 3).
    """
    # Get the bounding box of the alpha shape
    bbox_min, bbox_max = alpha_shape_mesh.bounds

    # Generate random points within the bounding box
    sampled_points = []
    while len(sampled_points) < num_samples:
        # Generate random points within the bounding box
        points = np.random.uniform(low=bbox_min, high=bbox_max, size=(num_samples * 2, 3) )
        inside = alpha_shape_mesh.contains(points)
        sampled_points.extend(points[inside])
        print(len(sampled_points))
    # Return the requested number of samples
    return np.array(sampled_points)


import random


def save_alpha_shape_as_ply(alpha_shape, file_path):
    """
    Save a Trimesh alpha shape as a PLY file.

    Args:
        alpha_shape: A Trimesh object representing the alpha shape.
        file_path: The output file path for the PLY file.
    """
    if not isinstance(alpha_shape, trimesh.Trimesh):
        raise TypeError("avislpha_shape must be a trimesh.Trimesh object.")

    # Save as PLYencoding='ascii'
    alpha_shape.export(file_path, file_type="ply", encoding="ascii")
    print(f"Alpha shape saved as PLY file at {file_path}")

def produce_alpha_contour(RCSB_ID, alpha):
    num_samples = 5000
    cifpath = AssetType.MMCIF.get_path(RCSB_ID)
    point_cloud = cif_to_point_cloud(cifpath)
    point_cloud = np.array(point_cloud)
    print("Produced point cloud for {}".format(RCSB_ID))

    print("Constructing alpha shape with a={}".format(alpha))
    alpha_shape = alphashape.alphashape(point_cloud, alpha)
    # alpha_shape.show()

    components = alpha_shape.split(only_watertight=False)  # Get all components
    alpha_shape_largest = max(components, key=lambda c: abs(c.volume))
    random.seed(10)

    print("Resampling {} points on initial alphashape M.".format(num_samples))
    new_points = sample_within_alpha_shape(alpha_shape_largest, num_samples)
    alpha_shape_renew = alphashape.alphashape(new_points, alpha)
    # alpha_shape_renew.show()

    path = AssetType.ALPHA_SHAPE.get_path(RCSB_ID)
    save_alpha_shape_as_ply(alpha_shape_renew, path)
    print("Saved to {}".format(path))
