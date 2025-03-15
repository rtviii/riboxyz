#!/usr/bin/env python3
import os
import json
import re
from collections import defaultdict

from ribctl.lib.libtax import Taxid

def determine_kingdom(structure_data):
    """
    Determine the kingdom of a ribosome structure.
    First check if it's mitochondrial, then use taxonomy ID to determine kingdom.
    """
    # Check if it's mitochondrial first
    if structure_data.get('mitochondrial', False):
        return "mitochondria"
    
    # Get the source organism IDs
    src_organism_ids = structure_data.get('src_organism_ids', [])
    if not src_organism_ids:
        return "unknown"
    
    # Use the first organism ID for classification
    taxid = src_organism_ids[0]
    
    try:
        
        # Constants for taxonomy IDs
        TAXID_BACTERIA = 2
        TAXID_EUKARYOTA = 2759
        TAXID_ARCHAEA = 2157
        
        if Taxid.is_descendant_of(TAXID_BACTERIA, taxid):
            return "bacteria"
        elif Taxid.is_descendant_of(TAXID_EUKARYOTA, taxid):
            return "eukarya"
        elif Taxid.is_descendant_of(TAXID_ARCHAEA, taxid):
            return "archaea"
        else:
            return "unknown"
    except Exception as e:
        print(f"Error determining kingdom for taxid {taxid}: {e}")
        return "unknown"

def get_species_name(structure_data):
    """Get the species name from the structure data."""
    if 'src_organism_names' in structure_data and structure_data['src_organism_names']:
        return structure_data['src_organism_names'][0]
    return "Unknown species"

def generate_kingdom_report():
    """Generate a report on the distribution of kingdoms among the mesh files."""
    mesh_dir = "wenjun_data_meshes"
    profile_dir = "wenjun_data_profiles"
    
    # Get all mesh files
    mesh_files = [f for f in os.listdir(mesh_dir) if f.endswith('.ply')]
    
    # Dictionaries to store results
    kingdom_counts = defaultdict(int)
    kingdom_ids = defaultdict(list)
    species_info = {}  # Map PDB IDs to species names
    
    for mesh_file in mesh_files:
        # Extract PDB ID from mesh filename
        match = re.match(r'([A-Za-z0-9]+)_NPET_MESH\.ply', mesh_file)
        if not match:
            print(f"Could not extract PDB ID from {mesh_file}")
            continue
        
        pdb_id = match.group(1)
        profile_file = os.path.join(profile_dir, f"{pdb_id}.json")
        
        # Check if corresponding profile exists
        if not os.path.exists(profile_file):
            print(f"Profile file not found for {pdb_id}")
            continue
        
        # Load profile data
        try:
            with open(profile_file, 'r') as f:
                structure_data = json.load(f)
            
            # Determine kingdom
            kingdom = determine_kingdom(structure_data)
            
            # Get species name
            species = get_species_name(structure_data)
            species_info[pdb_id] = species
            
            # Update counts and lists
            kingdom_counts[kingdom] += 1
            kingdom_ids[kingdom].append(pdb_id)
            
        except Exception as e:
            print(f"Error processing {pdb_id}: {e}")
    
    # Generate report
    print("\n--- Kingdom Distribution Report ---")
    print(f"Total structures analyzed: {sum(kingdom_counts.values())}")
    
    print("\nDistribution:")
    kingdoms = ["bacteria", "archaea", "eukarya", "mitochondria", "unknown"]
    for kingdom in kingdoms:
        if kingdom in kingdom_counts:
            count = kingdom_counts[kingdom]
            percentage = (count / sum(kingdom_counts.values())) * 100
            print(f"- {kingdom.capitalize()}: {count} structures ({percentage:.1f}%)")
    
    print("\nStructures by Kingdom:")
    for kingdom in kingdoms:
        if kingdom in kingdom_ids and kingdom_ids[kingdom]:
            print(f"\n{kingdom.capitalize()} ({len(kingdom_ids[kingdom])}):")
            
            # Group by species
            species_groups = defaultdict(list)
            for pdb_id in kingdom_ids[kingdom]:
                species = species_info.get(pdb_id, "Unknown species")
                species_groups[species].append(pdb_id)
            
            # Print species-wise grouping
            for species, ids in sorted(species_groups.items()):
                print(f"  {species} ({len(ids)}):")
                # Print IDs in groups of 5 for readability
                for i in range(0, len(ids), 5):
                    print("    " + ", ".join(sorted(ids[i:i+5])))

if __name__ == "__main__":
    generate_kingdom_report()
