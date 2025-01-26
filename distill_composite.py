import json
from collections import defaultdict
import os
from typing import List, Dict, Any, Tuple

from ribctl import ASSETS_PATH

def process_binding_sites(items: List[Dict[str, Any]]) -> Dict[str, Dict]:
    """
    Process ligand items and group them by ribosome_ligand_categories.
    
    Args:
        items: List of ligand dictionaries containing the data
        
    Returns:
        Dictionary with categories as keys and processed data as values
    """
    # Initialize result dictionary using defaultdict
    result = defaultdict(lambda: {"items": [], "composite_bsite": set()})
    
    for item in items:
        # Extract the relevant fields
        processed_item = {
            "chemicalId": item.get("chemicalId"),
            "chemicalName": item.get("chemicalName"),
            "purported_7K00_binding_site": item.get("purported_7K00_binding_site", [])
        }
        
        # Get categories, default to ["Other"] if none exists
        categories = item.get("ribosome_ligand_categories", ["Other"])
        if not categories:  # Handle empty list case
            categories = ["Other"]
            
        # Add item to each of its categories
        for category in categories:
            result[category]["items"].append(processed_item)
            
            # Update composite binding site
            if processed_item["purported_7K00_binding_site"]:
                # Convert binding site tuples to a consistent format for set operations
                binding_sites = set(tuple(site) for site in processed_item["purported_7K00_binding_site"])
                result[category]["composite_bsite"].update(binding_sites)
    
    # Convert sets back to sorted lists for JSON serialization
    final_result = {}
    for category, data in result.items():
        final_result[category] = {
            "items": data["items"],
            "composite_bsite": sorted(list(map(list, data["composite_bsite"])), 
                                   key=lambda x: (x[0], x[1]))
        }
    
    return final_result

def main():
    merged_ligand_data_path = os.path.join(ASSETS_PATH, 'ligands','merged_ligand_data.json')
    with open(merged_ligand_data_path, 'r') as f:
        data = json.load(f)
    
    # Handle both single item and list of items
    items = [data] if not isinstance(data, list) else data
    
    # Process the data
    result = process_binding_sites(items)
    
    # Write the result to a new JSON file
    with open('processed_ligands.json', 'w') as f:
        json.dump(result, f, indent=2)
    
    # Print some statistics
    print("\nProcessing completed!")
    print(f"Number of categories found: {len(result)}")
    for category, data in result.items():
        print(f"\n{category}:")
        print(f"  Number of items: {len(data['items'])}")
        print(f"  Size of composite binding site: {len(data['composite_bsite'])}")

if __name__ == "__main__":
    main()