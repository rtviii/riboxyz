import json
import os
from ribctl import ASSETS_PATH


def merge_ligand_data(binding_sites_dict, classyfire_list, chemical_details_list):
    """ Merge ligand data from three different sources based on matching identifiers. """
    merged_data = {}

    for classyfire_obj in classyfire_list:
        print(classyfire_obj)
        if "identifier" not in classyfire_obj:
            continue
        ligand_id = classyfire_obj["identifier"]
        merged_data[ligand_id] = classyfire_obj.copy()

    for chem_details in chemical_details_list:
        chemical_id = chem_details["chemicalId"]
        if chemical_id in merged_data:
            merged_data[chemical_id].update(chem_details)
        else:
            merged_data[chemical_id] = chem_details.copy()

    for ligand_id, binding_site in binding_sites_dict.items():
        if ligand_id in merged_data:
            merged_data[ligand_id]["purported_7K00_binding_site"] = binding_site
        else:
            merged_data[ligand_id] = {
                "identifier": ligand_id,
                "purported_7K00_binding_site": binding_site,
            }

    return list(merged_data.values())


if __name__ == "__main__":
    # Load real data from JSON files
    binding_sites_path = os.path.join(
        ASSETS_PATH, "ligands", "7K00_ligand_predictions.json"
    )
    classyfire_path = os.path.join(ASSETS_PATH, "ligands", "ligand_classes.json")
    chemical_details_path = os.path.join(ASSETS_PATH, "ligands", "ribxz_ligands.json")

    # Open files with utf-8-sig encoding to handle BOM
    with open(binding_sites_path, "r", encoding="utf-8-sig") as f:
        binding_sites = json.load(f)

    with open(classyfire_path, "r", encoding="utf-8-sig") as f:
        classyfire_data = json.load(f)

    with open(chemical_details_path, "r", encoding="utf-8-sig") as f:
        chemical_details = json.load(f)

    # Merge the data
    merged_results = merge_ligand_data(binding_sites, classyfire_data, chemical_details)

    # Save merged results
    output_path = os.path.join(ASSETS_PATH, "ligands", "merged_ligand_data.json")
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(merged_results, f, indent=2)

    print(f"Merged data saved to {output_path}")
