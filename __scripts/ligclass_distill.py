import json
import os
from typing import Dict, List, Set
from ribctl import ASSETS_PATH

# Primary categories we care about
PRIMARY_CATEGORIES = [
    "Aminoglycosides",
    "Macrolides",
    "Tetracyclines",
    "Oxazolidinones"
]

def normalize_string(s: str) -> str:
    """Normalize string for comparison"""
    return s.lower().strip()

def extract_categories(ligand: dict) -> List[str]:
    """Extract relevant categories from ligand data"""
    if not ligand:
        return []
    
    categories: Set[str] = set()
    
    def process_category(category):
        if not category:
            return
        
        normalized_category = normalize_string(category)
        primary_match = next(
            (cat for cat in PRIMARY_CATEGORIES 
             if normalize_string(cat) in normalized_category
             or normalized_category in normalize_string(cat)),
            None
        )
        
        if primary_match:
            categories.add(primary_match)
    
    # Process all category sources
    if ligand.get('direct_parent'):
        process_category(ligand['direct_parent'].get('name'))
    
    for parent in ligand.get('alternative_parents', []):
        process_category(parent.get('name'))
    
    for ancestor in ligand.get('ancestors', []):
        process_category(ancestor)
    
    for node in ligand.get('intermediate_nodes', []):
        process_category(node.get('name'))
    
    for descriptor in ligand.get('external_descriptors', []):
        for annotation in descriptor.get('annotations', []):
            process_category(annotation)
    
    return list(categories)


def has_sufficient_residues(ligand: dict, min_residues: int = 10) -> bool:
    """
    Check if ligand has sufficient number of residues in binding site.
    
    Args:
        ligand (dict): Ligand data dictionary
        min_residues (int): Minimum number of residues required (default: 10)
    
    Returns:
        bool: True if ligand has sufficient residues, False otherwise
    """
    binding_sites = ligand.get('purported_7K00_binding_site', [])
    return len(binding_sites) >= min_residues

def create_simplified_data():
    """Create simplified version of merged data"""
    # Load merged data
    merged_data_path = os.path.join(ASSETS_PATH, 'ligands', 'merged_ligand_data.json')
    with open(merged_data_path, 'r', encoding='utf-8') as f:
        merged_data = json.load(f)
    
    # Create simplified version
    simplified_data = []
    
    for ligand in merged_data:
        # Skip ligands with insufficient residues
        if not has_sufficient_residues(ligand):
            continue
            
        # Extract categories
        tags = extract_categories(ligand)
        
        # Create simplified ligand object
        simple_ligand = {
            "chemicalId"                 : ligand.get('chemicalId') or ligand.get('identifier'),
            "chemicalName"               : ligand.get('chemicalName') or ligand.get('pdbx_description'),
            "purported_7K00_binding_site": ligand.get('purported_7K00_binding_site'),
            "tags"                       : tags
        }
        
        simplified_data.append(simple_ligand)

    return merged_data_path

if __name__ == "__main__":
    import sys
    import tempfile
    import shutil

    # Get input file path from command line or use default
    input_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ASSETS_PATH, 'ligands', 'ligands_data_7k00_demo.json')
    
    # Load existing data
    with open(input_path, 'r', encoding='utf-8') as f:
        data = json.load(f)
    
    # Create simplified data
    simplified_data = []
    for ligand in data:
        # Skip ligands with insufficient residues
        if not ligand.get('purported_7K00_binding_site') or not has_sufficient_residues(ligand) :
            continue
            
        # Extract categories
        tags = extract_categories(ligand)
        
        # Create simplified ligand object
        simple_ligand = {
            "chemicalId": ligand.get('chemicalId') or ligand.get('identifier'),
            "chemicalName": ligand.get('chemicalName') or ligand.get('pdbx_description'),
            "purported_7K00_binding_site": ligand.get('purported_7K00_binding_site'),
            "tags": tags
        }
        
        simplified_data.append(simple_ligand)
    
    # Write to temporary file first (safer than direct overwrite)
    with tempfile.NamedTemporaryFile(mode='w', delete=False, encoding='utf-8') as tmp:
        json.dump(simplified_data, tmp, indent=2)
    
    # Replace original file with new file
    shutil.move(tmp.name, input_path)
    
    print(f"File processed and updated in place: {input_path}")
    print(f"Filtered from {len(data)} to {len(simplified_data)} ligands")