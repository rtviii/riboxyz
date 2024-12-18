import argparse
from pprint import pprint
import sys
sys.path.append('/home/rtviii/dev/riboxyz')
from matplotlib import pyplot as plt
import numpy as np
from sklearn.cluster import DBSCAN
from mesh_generation.mesh_bbox_extraction import ( encode_atoms, open_tunnel_csv, parse_struct_via_bbox, parse_struct_via_centerline)
from mesh_paths import *
from mesh_generation.mesh_full_pipeline import pipeline
from mesh_generation.mes_visualization import plot_with_landmarks


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


        
def main():

    # ? ---------- Params ------------
    parser = argparse.ArgumentParser()
    # Add command-line arguments

    parser.add_argument( "--full_pipeline",   action='store_true')

    # visualization options
    parser.add_argument( "--result",   action='store_true')
    parser.add_argument( "--gif_intermediates",   action='store_true')


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

    #! Interactive / processing
    parser.add_argument( "--trim",  action='store_true',required=False)
    parser.add_argument( "--nascent_chain",  action='store_true',required=False)
    parser.add_argument( "--cluster_manual",  action='store_true',required=False)

    args          = parser.parse_args()
    RCSB_ID       = args.rcsb_id.upper()

    if args.full_pipeline:
        pipeline(RCSB_ID, args)
        exit(0)
    if args.result:
        plot_with_landmarks(RCSB_ID,None,args.gif_intermediates, "{}.result.gif".format(RCSB_ID))



if __name__ == "__main__":
    main()


