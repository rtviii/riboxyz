import gzip
import os
import numpy as np
import requests
from scipy.spatial.distance import cdist


def midpoint(p1, p2):
    return (p1 + p2) / 2

def find_closest_pair_two_sets(points1, points2):
    """
    Find the pair of points (one from each set) that are closest to each other.
    
    Parameters:
    points1: numpy array of shape (N, 3) for first set of 3D points
    points2: numpy array of shape (M, 3) for second set of 3D points
    
    Returns:
    tuple: (point1, point2, distance) where point1 and point2 are the closest points
           and distance is their Euclidean distance
    """
    # Convert inputs to numpy arrays if they aren't already
    points1 = np.asarray(points1)
    points2 = np.asarray(points2)
    
    # Calculate pairwise distances between all points
    distances = cdist(points1, points2)
    
    # Find the indices of the minimum distance
    i, j = np.unravel_index(distances.argmin(), distances.shape)
    
    # Get the closest points and their distance
    closest_point1 = points1[i]
    closest_point2 = points2[j]
    min_distance = distances[i, j]
    
    return closest_point1, closest_point2

async def download_unpack_place(struct_id: str) -> None:
    BASE_URL = "http://files.rcsb.org/download/"
    FORMAT   = ".cif.gz"

    structid     = struct_id.upper()
    url          = BASE_URL + structid + FORMAT
    compressed   = requests.get(url).content
    decompressed = gzip.decompress(compressed)

    structfile = os.path.join(
        os.environ["RIBETL_DATA"],
        structid,
        structid + ".cif")

    with open(structfile, "wb") as f:
        f.write(decompressed)

