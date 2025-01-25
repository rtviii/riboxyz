import json
import os

from ribctl import ASSETS_PATH

# Load the JSON data
with open(os.path.join(ASSETS_PATH, 'ligands','merged_ligand_data.json'), 'r') as file:
    ligands = json.load(file)
    print(ligands)


# Define the classification rules
def classify_ligand(ligand):
    # Safeguard: If ligand is None, return an empty list of classes
    if ligand is None:
        return []
    
    classes = []
    
    # Safely flatten the relevant fields into an array
    flattened_fields = [
        (ligand.get('kingdom') or {}).get('name', ''),  # Safeguard for 'kingdom'
        (ligand.get('superclass') or {}).get('name', ''),  # Safeguard for 'superclass'
        (ligand.get('class') or {}).get('name', ''),  # Safeguard for 'class'
        (ligand.get('subclass') or {}).get('name', ''),  # Safeguard for 'subclass'
        (ligand.get('direct_parent') or {}).get('name', '')  # Safeguard for 'direct_parent'
    ]
    
    # Add intermediate nodes (if present)
    flattened_fields.extend([(node or {}).get('name', '') for node in ligand.get('intermediate_nodes', [])])
    
    # Add alternative parents (if present)
    flattened_fields.extend([(parent or {}).get('name', '') for parent in ligand.get('alternative_parents', [])])
    
    # Add predicted ChEBI terms (if present)
    flattened_fields.extend([term.split(' (')[0] for term in ligand.get('predicted_chebi_terms', [])])
    
    # Add external descriptors (if present)
    for descriptor in ligand.get('external_descriptors', []):
        flattened_fields.extend(descriptor.get('annotations', []))
    
    # Define the categories to search for
    categories = [
        'Aminoglycosides', 'Macrolides', 'Tetracyclines', 'Oxazolidinones', 'Ketolides','Aminocyclitols',
        'Phenicols',  'Lincosamides',
         'Pleuromutilins',   'Streptogramins'
    ]
    
    # Check for each category in the flattened fields
    for category in categories:
        if any(category.lower() in field.lower() for field in flattened_fields if field):  # Skip empty fields
            classes.append(category)
    
    return classes

# Add categories to each ligand and write to a new file
for ligand in ligands:
    if ligand is not None:  # Skip None entries
        categories = classify_ligand(ligand)
        ligand['ribosome_ligand_categories'] = categories

# Write the updated data to a new JSON file
output_file = 'merged_ligand_data_with_categories.json'
with open(output_file, 'w') as file:
    json.dump(ligands, file, indent=4)

print(f"Updated data has been written to {output_file}")