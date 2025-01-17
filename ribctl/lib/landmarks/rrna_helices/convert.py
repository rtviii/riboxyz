import csv
from pprint import pprint
import sys
import json

sys.path.append(
    "/home/rtviii/dev/riboxyz"
)  #! hack until ribctl is a separate pypi project (after that just pip install ribctl)
from io import StringIO

from ribctl.ribosome_ops import RibosomeOps


def convert_csv_to_json(csv_text):
    # Initialize the result dictionary
    result = {}

    # Read CSV
    reader = csv.DictReader(StringIO(csv_text))

    for row in reader:
        # Extract chain (23S, 16S, 5S)
        chain = row["ResNum"].split(":")[0]

        # Extract helix name
        helix_name = row["HelixName"]

        # Extract first range only (everything before the semicolon)
        first_range = row["ResNum"].split(";")[0]

        # Extract numbers from the range using string operations
        start, end = map(int, first_range.split("(")[1].split(")")[0].split("-"))

        # Initialize chain if it doesn't exist
        if chain not in result:
            result[chain] = {}

        # Add helix with range
        result[chain][helix_name] = [[start, end]]

    return result


# You can either read from file:
# name = 'thermus_1vy4'
# # name = 'yeast_4v88'
# with open(f'{name}.csv', 'r') as f:
#     csv_text = f.read()
#     result = convert_csv_to_json(csv_text)
#     with open(f'{name}.json', 'w') as f:
#         f.write(json.dumps(result, indent=2))
#     print(json.dumps(result, indent=2))

with open(
    "/home/rtviii/dev/riboxyz/ribctl/lib/landmarks/rrna_helices/ecoli_7K00.json", "r"
) as f:
    reference_helices = json.load(f)

rcsb_id = "7K00"
profile = RibosomeOps(rcsb_id).profile
p = profile.get_nomenclature_map()
reversed_nom = {_[1][0]: _[0] for _ in p.items() if len(_[1])}

aaid_16s = reversed_nom["rRNA_16S"]
aaid_23s = reversed_nom["rRNA_23S"]
aaid_5s  = reversed_nom["rRNA_5S"]

print(aaid_16s, aaid_23s, aaid_5s)


new = {}
for chain in ["16S", "23S", "5S"]:
    if chain == "16S":
        aaid = aaid_16s
        new[aaid] = {"polymer_class": "16SrRNA", "helices": {}}
    elif chain == "23S":
        aaid = aaid_23s
        new[aaid] = {"polymer_class": "23SrRNA", "helices": {}}
    elif chain == "5S":
        aaid = aaid_5s
        new[aaid] = {"polymer_class": "5SrRNA", "helices": {}}
    else:
        raise ValueError("Invalid chain")

    for helix_name, ranges in reference_helices[chain].items():
        new[aaid]["helices"][helix_name] = [*ranges[0]]
pprint(new)

s23 = RibosomeOps(rcsb_id).get_LSU_rRNA()
pprint(s23.entity_poly_seq_one_letter_code.__len__())
pprint(s23.entity_poly_seq_one_letter_code_can.__len__())

with open(f"{rcsb_id}_rrna_helices.json", "w") as f:
    f.write(json.dumps(new, indent=2))