#!/usr/bin/env python
import os
import sys
import glob
import json
import argparse
import numpy as np
from pathlib import Path
from tqdm import tqdm
import scipy.spatial.distance as distance

from ribctl.ribosome_ops import RibosomeOps
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.schema.types_ribosome import ConstrictionSite, PTCInfo

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

def process_structure(rcsb_id):
    """Process a single structure and extract landmark data"""
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
        
        # Load PTC information
        ptc_path = AssetType.PTC.get_path(rcsb_id)
        if ptc_path.exists():
            with open(ptc_path, 'r') as f:
                ptc_info = PTCInfo.model_validate(json.load(f))
            landmarks["ptc"] = ptc_info.location
        else:
            print(f"  Warning: No PTC data found for {rcsb_id}")
        
        # Load constriction site information
        constriction_path = AssetType.CONSTRICTION_SITE.get_path(rcsb_id)
        if constriction_path.exists():
            with open(constriction_path, 'r') as f:
                constriction_info = ConstrictionSite.model_validate(json.load(f))
            landmarks["constriction"] = constriction_info.location
        else:
            print(f"  Warning: No constriction site data found for {rcsb_id}")
            # If no constriction site, we can't find closest residues
            return landmarks
        
        # Find uL4 and uL22 protein chains
        profile = ro.profile
        ul4_chain = None
        ul22_chain = None
        
        # Search through proteins to find uL4 and uL22
        for protein in profile.proteins:
            for nomenclature in protein.nomenclature:
                if nomenclature.value == "uL4":
                    ul4_chain = protein.auth_asym_id
                elif nomenclature.value == "uL22":
                    ul22_chain = protein.auth_asym_id
        
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
    
    return landmarks

def main():

    parser = argparse.ArgumentParser(description='Extract landmark data from ribosome structures')
    parser.add_argument('--data-dir', type=str, help='Directory containing structure data')
    parser.add_argument('--output-dir', type=str, help='Directory to save output files')
    parser.add_argument('--ids', type=str, nargs='+', help='Specific structure IDs to process')
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
    
    # Process each structure
    for rcsb_id in tqdm(structure_ids):
        output_file = os.path.join(output_dir, f"{rcsb_id}_landmarks_pts.json")
        
        # Skip if file already exists (unless specific IDs were requested)
        if os.path.exists(output_file) and not args.ids:
            print(f"Skipping {rcsb_id}, file already exists")
            continue
        
        # Process the structure
        landmarks = process_structure(rcsb_id)
        
        # Save to file
        with open(output_file, 'w') as f:
            json.dump(landmarks, f, indent=2)
        
        print(f"Saved landmark data for {rcsb_id}")
    
    print("Done!")

if __name__ == "__main__":
    main()