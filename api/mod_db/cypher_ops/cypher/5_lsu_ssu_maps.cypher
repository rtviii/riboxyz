call apoc.load.json('file:///import/riboxyz_seed_data/ban-pfam-map-ssu.json') yield value
unwind(keys(value)) as ssuclass merge (ncs:NomenclatureClass {class_id:ssuclass});

call apoc.load.json('file:///import/riboxyz_seed_data/ban-pfam-map-lsu.json') yield value
unwind(keys(value)) as  lsuclass merge (ncl:NomenclatureClass {class_id:lsuclass});