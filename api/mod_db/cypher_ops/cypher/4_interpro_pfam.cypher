WITH [1,2,3,4] AS coll
UNWIND coll AS x WITH x
CALL apoc.load.json(apoc.text.format("file:///import/riboxyz_seed_data/pfam-interpro-%s.json", [x])  ) yield value as entry
with entry.metadata as datum
with datum where datum.integrated is not null
merge (inode:InterProFamily{family_id: datum.integrated})
merge (pnode:PFAMFamily{family_id: datum.accession, family_type:datum.type})
merge (inode)-[:mp_InterPro_PFAM]-(pnode)