WITH [1,2,3,4] AS coll
UNWIND coll AS x
WITH x
CALL apoc.load.json(apoc.text.format('file:///import/riboxyz_seed_data/interpro-go-%s.json', [x])) yield value as go
merge (inode:InterProFamily{family_id:go.InterPro})
merge (gonode:GOClass{go_class:go.GO})
on create set gonode.annotation = go.GO_annotation
merge (inode)-[:mp_InterPro_GO{annotation:go.interpro_class}]-(gonode)
