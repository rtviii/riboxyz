
Get ligandcounts.
```
find . -type f -name '[A-Z0-9][A-Z0-9][A-Z0-9][A-Z0-9].json' -exec sh -c 'jq ".nonpolymeric_ligands[] | {chemicalId, number_of_instances}" "{}" >> ligandcounts.txt' \;
```
