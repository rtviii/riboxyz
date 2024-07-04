from pprint import pprint
import sys
sys.path.append("/home/rtviii/dev/riboxyz/ribctl")
from ribctl.etl.etl_assets_ops import RibosomeOps





dataset = [
    ["5NO2", "4V61"],  #<
    ["3J9M", "4CE4"],  #<
    ["5ND9", "6FXC"],  #<
    ["5JUU", "5AJ4"],
    ["5X8P", "5MRF"],
    ["5ADY", "5KCS"],  #<
    ["3J7Y", "5MLC"],
    ["5XXB", "5MLC"],
    ["3J7Y", "3J78"],  #< somehow
    # ["5FKU", "5M1S"],
]

ix = {}
for j,i in enumerate( dataset ):
    print("PAIR # {} | {}".format(j+1, i))
    n1 = [ ]
    for k in [kvp[1]['nomenclature'] for kvp in  RibosomeOps(i[0]).nomenclature_table().items() ]:
        n1 = [*n1, *k]
    n2 = [ ]
    for s in [ kvp[1]['nomenclature'] for kvp in  RibosomeOps(i[1]).nomenclature_table().items() ]:
        n2= [*n2, *s]

    if j==0:
        print(sorted(n1))
        print(sorted(n2))
    ix.update( { f"pair_{j+1}_{i[0]}_{i[1]}": set.intersection(set(n1), set(n2))  })

# pprint(ix)


for c in ix["pair_1_5NO2_4V61"] :
    print(RibosomeOps("4v61").get_poly_by_polyclass(c).auth_asym_id)