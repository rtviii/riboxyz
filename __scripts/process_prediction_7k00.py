import os
import json
from typing import Dict, List, Tuple
from pydantic import BaseModel
from pathlib import Path

from ribctl.lib.schema.types_binding_site import LigandTransposition

def process_prediction_files(directory: str) -> Dict[str, Dict[str, List[Tuple[int, str]]]]:
    """
    Process prediction JSON files and organize by ligand and chain.
    Returns: {ligand_id: {auth_asym_id: [(auth_seq_id, auth_asym_id), ...]}}
    """
    result = {}
    
    for filename in os.listdir(directory):
        if not filename.endswith('.json') or not filename.startswith('7K00_LIG'):
            continue
            
        filepath = os.path.join(directory, filename)
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
                # Validate data structure
                transposition = LigandTransposition.model_validate(data)
                
                ligand_id = transposition.purported_binding_site.ligand
                if ligand_id not in result:
                    result[ligand_id] = {}
                    
                # Process each chain's predictions
                for chain in transposition.purported_binding_site.chains:
                    chain_id = chain.auth_asym_id
                    if chain_id not in result[ligand_id]:
                        result[ligand_id][chain_id] = []
                        
                    # Add bound residues
                    for residue in chain.bound_residues:
                        result[ligand_id][chain_id].append(
                            (residue.auth_seq_id, residue.auth_asym_id)
                        )
                        
        except Exception as e:
            print(f"Error processing {filename}: {str(e)}")
            continue
            
    return result

if __name__ == "__main__":
    directory = "/home/rtviii/dev/RIBETL_DATA/7K00"
    predictions = process_prediction_files(directory)
    
    # Save organized results
    output_path = "organized_predictions.json"
    with open(output_path, 'w') as f:
        json.dump(predictions, f, indent=2)
    
    print(f"Processed predictions saved to {output_path}")