# How to use Classifyre
Basically i have a script that grab all my ligands and their "SMILES" codes, collates them into a single parametrized url request and hits the  http://classyfire.wishartlab.com/queries API. Param template below. The query result returns to you the id of your "ticket". Save it.

What you have to do later (after ~10min-1h) is go to their website with `http://classyfire.wishartlab.com/queries/YOUR_ID` and look up the results. Some will fail, but the majority should be classified. I emailed them to see if the source is avaialble but not holding my breath. Useful tool though.



```python
import sys
sys.path.append("/home/rtviii/dev/riboxyz")
from pprint import pprint
from rdkit import Chem
from rdkit.Chem import AllChem
from neo4j_ribosome.db_lib_reader import Neo4jReader
from ribctl.lib.schema.types_ribosome import NonpolymericLigand

def collect_all_ligands_with_smiles():
    reader      = Neo4jReader()
    ligs        = reader.list_ligands(nodes_only=True)[0][0]
    print("Got ligands:", len(ligs))
    ligands     = list(map(NonpolymericLigand.model_validate, ligs))
    query_input = ""
    lig:NonpolymericLigand
    for lig in ligands:
        query_input = query_input + f"{lig.chemicalId}\t{lig.SMILES_stereo}\n"
        # query_input = query_input + f"{lig.chemicalId}\t {lig.chemicalName}\n"
    return query_input
print(collect_all_ligands_with_smiles())


def init_classification_class():
  import requests
  url = "http://classyfire.wishartlab.com/queries"

  headers = {
      "Accept": "application/json",
      "Content-Type": "application/json"
  }

  data = {
      "label"      : "All ligands | SMiles stereo",
      # "query_input": "MOL1\\tCCCOCC\\nMOL2\\tCOCC=CCCC",
      "query_input": collect_all_ligands_with_smiles(),
      "query_type" : "STRUCTURE"
  }

  response = requests.post(url, json=data, headers=headers)

  print(response.status_code)
  print(response.json())


# init_classification_class()
# Example for http://classyfire.wishartlab.com/queries
# curl  -H "Accept: application/json" -H "Content-type: application/json" -d '{"label":"test", "query_input":"PAR\tCC1[C@H]([C@@H]([C@H]([C@@H]([C@H]1N)O[C@@H]2[C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)N)O[C@H]3[C@@H]([C@@H]([C@H](O3)CO)O[C@@H]4[C@@H]([C@H]([C@@H]([C@@H](O4)CN)O)O)N)O)O)N\nKSG\tC[C@@H]1[C@H](C[C@@H]([C@H](O1)OC2[C@@H]([C@H](C([C@@H]([C@@H]2O)O)O)O)O)N)N=C(C(=O)O)N", "query_type":"STRUCTURE"}' -X POST http://classyfire.wishartlab.com/queries
12123363
```