


#!/usr/bin/env python
import os
import sys
import glob
import json
import argparse
import numpy as np
from pathlib import Path
import multiprocessing as mp
import scipy.spatial.distance as distance
from functools import partial

from Bio.PDB.Chain import Chain
from ribctl.lib.types.polymer.base import CytosolicProteinClass, MitochondrialProteinClass
from ribctl.lib.utils import find_closest_pair_two_sets, midpoint
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.landmarks.ptc_via_trna import PTC_location


def get_constriction(rcsb_id: str) -> np.ndarray:
    """Calculate constriction site based on uL4 and uL22 protein chains"""
    ro = RibosomeOps(rcsb_id)
    is_mitochondrial = ro.profile.mitochondrial
    if is_mitochondrial:
        uL4 = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL4m)
        uL22 = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL22m)
    else:
        uL4 = ro.get_poly_by_polyclass(CytosolicProteinClass.uL4)
        uL22 = ro.get_poly_by_polyclass(CytosolicProteinClass.uL22)
    if uL4 is None or uL22 is None:
        raise ValueError(f"Could not find uL4 or uL22 in {rcsb_id}")
    structure = ro.assets.biopython_structure()
    
    uL4_c: Chain = structure[0][uL4.auth_asym_id]
    uL22_c: Chain = structure[0][uL22.auth_asym_id]
    uL4_coords = [r.center_of_mass() for r in uL4_c.child_list]
    uL22_coords = [r_.center_of_mass() for r_ in uL22_c.child_list]
    return midpoint(*find_closest_pair_two_sets(uL4_coords, uL22_coords))


def get_closest_residues(residues, target_location, n=15):
    """Get the n residues with centers closest to the target location"""
    residue_centers = []
    
    for residue in residues:
        try:
            # Calculate center of residue
            coords = np.array([atom.coord for atom in residue.get_atoms()])
            if len(coords) > 0:
                center = coords.mean(axis=0)
                residue_centers.append(center)
        except Exception:
            continue
    
    if not residue_centers:
        return np.array([])
    
    # Calculate distances to target location
    residue_centers = np.array(residue_centers)
    distances = distance.cdist([target_location], residue_centers)[0]
    
    # Get indices of n closest residues
    closest_indices = np.argsort(distances)[:n]
    
    # Return the coordinates of the centers of the n closest residues
    return residue_centers[closest_indices]


def process_structure(rcsb_id, output_dir=None, skip_existing=True):
    """Process a single structure and extract landmark data with zero transformations"""
    # Check if output file already exists (if output_dir and skip_existing are provided)
    if output_dir and skip_existing:
        output_file = os.path.join(output_dir, f"{rcsb_id}_landmarks_pts.json")
        if os.path.exists(output_file):
            print(f"Skipping {rcsb_id}, file already exists")
            return None
    
    print(f"Processing {rcsb_id}...")
    
    landmarks = {
        "ptc": None,
        "constriction": None,
        "uL4": None,
        "uL22": None
    }
    
    try:
        # Initialize RibosomeOps
        ro = RibosomeOps(rcsb_id)
        
        # Compute PTC location using the function from ptc_via_trna
        try:
            ptc_info = PTC_location(rcsb_id)
            landmarks["ptc"] = ptc_info.location
            print(f"  Found PTC at {landmarks['ptc']}")
        except Exception as e:
            print(f"  Error finding PTC: {e}")
            return rcsb_id, landmarks
        
        # Compute constriction site location using the provided function
        try:
            constriction_location = get_constriction(rcsb_id)
            landmarks["constriction"] = constriction_location.tolist()
            print(f"  Found constriction site at {landmarks['constriction']}")
        except Exception as e:
            print(f"  Error finding constriction site: {e}")
            return rcsb_id, landmarks
        
        # Find uL4 and uL22 protein chains
        profile = ro.profile
        ul4_chain = None
        ul22_chain = None
        

        if ro.profile.mitochondrial:
            ul4_chain = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL4m).auth_asym_id
            ul22_chain = ro.get_poly_by_polyclass(MitochondrialProteinClass.uL22m).auth_asym_id
        else:
            ul4_chain = ro.get_poly_by_polyclass(CytosolicProteinClass.uL4).auth_asym_id
            ul22_chain = ro.get_poly_by_polyclass(CytosolicProteinClass.uL22).auth_asym_id

        # Process protein chains if found
        model = ro.assets.biopython_structure()[0]
        constriction_location = np.array(landmarks["constriction"])
        
        # Process uL4 chain
        if ul4_chain and ul4_chain in model.child_dict:
            ul4_residues = list(model[ul4_chain].get_residues())
            ul4_points = get_closest_residues(ul4_residues, constriction_location, 15)
            
            if len(ul4_points) > 0:
                landmarks["uL4"] = ul4_points.tolist()
            else:
                print(f"  Warning: No valid uL4 residues found for {rcsb_id}")
        else:
            print(f"  Warning: No uL4 chain found for {rcsb_id}")
        
        # Process uL22 chain
        if ul22_chain and ul22_chain in model.child_dict:
            ul22_residues = list(model[ul22_chain].get_residues())
            ul22_points = get_closest_residues(ul22_residues, constriction_location, 15)
            
            if len(ul22_points) > 0:
                landmarks["uL22"] = ul22_points.tolist()
            else:
                print(f"  Warning: No valid uL22 residues found for {rcsb_id}")
        else:
            print(f"  Warning: No uL22 chain found for {rcsb_id}")
    
    except Exception as e:
        print(f"  Error processing {rcsb_id}: {str(e)}")
    
    # Save to file if output_dir is provided
    if output_dir:
        output_file = os.path.join(output_dir, f"{rcsb_id}_landmarks_pts.json")
        with open(output_file, 'w') as f:
            json.dump(landmarks, f, indent=2)
        print(f"Saved landmark data for {rcsb_id}")
    
    return rcsb_id, landmarks


def process_worker(rcsb_id, output_dir, skip_existing):
    """Worker function for multiprocessing"""
    try:
        result = process_structure(rcsb_id, output_dir, skip_existing)
        if result:
            return result[0]  # Return rcsb_id if processed successfully
        return None  # Return None if skipped
    except Exception as e:
        print(f"Error in worker processing {rcsb_id}: {str(e)}")
        return None


def main():
    parser = argparse.ArgumentParser(description='Extract landmark data from ribosome structures')
    parser.add_argument('--data-dir', type=str, help='Directory containing structure data')
    parser.add_argument('--output-dir', type=str, help='Directory to save output files')
    parser.add_argument('--ids', type=str, nargs='+', help='Specific structure IDs to process')
    parser.add_argument('--processes', type=int, default=mp.cpu_count(), 
                      help='Number of processes to use (default: number of CPU cores)')
    parser.add_argument('--force', action='store_true', 
                      help='Force processing even if output file exists')
    args = parser.parse_args()
    
    # Set default data directory if not provided
    data_dir = args.data_dir if args.data_dir else os.environ.get("RIBETL_DATA", ".")
    
    # Set default output directory if not provided
    output_dir = args.output_dir if args.output_dir else os.path.join(data_dir, "landmarks")
    os.makedirs(output_dir, exist_ok=True)
    
    # Get structure IDs
    if args.ids:
        structure_ids = args.ids
    else:
        # Find all structure directories that have NPET meshes
        mesh_files = glob.glob(os.path.join(data_dir, "*/*_NPET_MESH.ply"))
        structure_ids = [os.path.basename(f).split('_')[0] for f in mesh_files]
        structure_ids = list(set(structure_ids))  # Remove duplicates
    
    print(f"Found {len(structure_ids)} structures to process")
    
    # Adjust number of processes if needed
    num_processes = min(args.processes, len(structure_ids))
    if num_processes < 1:
        num_processes = 1
    
    # Skip check based on force flag
    skip_existing = not args.force
    
    print(f"Using {num_processes} processes")
    
    # Process structures in parallel
    if num_processes > 1:
        # Create a partial function with fixed arguments
        worker_func = partial(process_worker, output_dir=output_dir, skip_existing=skip_existing)
        
        # Create pool and map worker function to structure IDs
        with mp.Pool(processes=num_processes) as pool:
            processed_ids = list(pool.map(worker_func, structure_ids))
        
        # Count successfully processed structures (non-None results)
        successful = [pid for pid in processed_ids if pid is not None]
        print(f"Successfully processed {len(successful)} structures")
    else:
        # Process sequentially if only one process requested
        for rcsb_id in structure_ids:
            process_structure(rcsb_id, output_dir, skip_existing)
    
    print("Done!")


if __name__ == "__main__":
    main()
