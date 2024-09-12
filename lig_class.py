from rdkit import Chem
from rdkit.Chem import AllChem
import requests
import json

def get_compound_class(smiles):
    # Generate InChI from SMILES using RDKit
    mol   = Chem.MolFromSmiles(smiles)
    inchi = Chem.MolToInchi(mol)

    # Use ClassyFire API to get classification
    url = "http://classyfire.wishartlab.com/queries"
    payload = {"query_input": inchi, "query_type": "INCHI"}
    headers = {"Content-Type": "application/json"}
    
    response = requests.post(url, data=json.dumps(payload), headers=headers)
    data     = response.json()
    
    # Get the query ID and wait for results
    query_id = data['id']
    result_url = f"http://classyfire.wishartlab.com/queries/{query_id}.json"
    
    while True:
        response = requests.get(result_url)
        data = response.json()
        if data['classification_status'] == 'Done':
            break
    
    classification = data['entities'][0]['direct_parent']
    return classification

# Example usage
smiles = "CC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)N)O[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O[C@@H]4[C@@H]([C@H]([C@@H]([C@@H](O4)CN)O)O)N)O)O)N"  # SMILES for Paromomycin
compound_class = get_compound_class(smiles)
print(f"The class of the compound is: {compound_class}")