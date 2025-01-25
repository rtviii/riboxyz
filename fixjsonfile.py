import os
import json
import re
from ribctl import ASSETS_PATH

def fix_smiles_json():
    chemical_details_path = os.path.join(ASSETS_PATH, 'ligands', 'ribxz_ligands.json')
    
    with open(chemical_details_path, 'r', encoding='utf-8-sig') as f:
        content = f.read()
    
    # Find all SMILES strings and fix their escapes
    pattern = r'"(SMILES|SMILES_stereo)"\s*:\s*"([^"]+)"'
    
    def fix_smiles(match):
        key = match.group(1)
        smiles = match.group(2)
        # Remove all backslashes from SMILES
        fixed_smiles = smiles.replace('\\', '')
        return f'"{key}": "{fixed_smiles}"'
    
    fixed_content = re.sub(pattern, fix_smiles, content)
    
    # Write fixed content
    fixed_path = os.path.join(ASSETS_PATH, 'ligands', 'ribxz_ligands_fixed.json')
    with open(fixed_path, 'w', encoding='utf-8') as f:
        f.write(fixed_content)
    
    # Test if it works
    try:
        with open(fixed_path, 'r', encoding='utf-8') as f:
            test_load = json.load(f)
            print("Successfully fixed and validated JSON!")
    except json.JSONDecodeError as e:
        print(f"Error: {e}")
        error_pos = int(str(e).split('char ')[-1].split(')')[0])
        print("Context around error:")
        print(fixed_content[error_pos-50:error_pos+50])
    
    return fixed_path

if __name__ == "__main__":
    fixed_path = fix_smiles_json()
    print(f"File saved to: {fixed_path}")
