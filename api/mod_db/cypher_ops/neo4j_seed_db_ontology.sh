#!/usr/bin/bash

# Make sure to set database.files.import variable in neo4j.conf correctly. that's what 'file:///' reads from by default
# (that's where file:// reads from by default : /var/lib/neo4j/import)
# the base data can be found at https://github.com/rtviii/neo4j-basedata.git

# You must set the database name in the script : "export RIBOXYZ_DB_NAME= XXX"

# Create constraints
./cypher_exec 1_constraints.cypher
./cypher_exec 2_interpro_base.cypher
./cypher_exec 3_interpro_go.cypher
./cypher_exec 4_interpro_pfam.cypher
./cypher_exec 5_lsu_ssu_maps.cypher
./cypher_exec 6_rna_classes.cypher


# # INTERPRO families
# echo 'call apoc.load.json("file:///ontology/data/interpro.json") yield value
# with value as v
# merge (q:InterProFamily{ family_id:KEYS(v)[0],type:v[KEYS(v)[0]].type,description:v[KEYS(v)[0]].name})' | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;


# # INTERPRO-GeneOntology connections
# for i in 1 2 3 4;
# do 
# 	echo "call apoc.load.json(\"file:///ontology/data/interpro-go/part$i.json\") yield value as go
# 	merge (inode:InterProFamily{family_id:go.InterPro})
# 	merge (gonode:GOClass{go_class:go.GO})
# 	on create set gonode.annotation = go.GO_annotation
# 	merge (inode)-[:interpro_go{annotation:go.interpro_class}]-(gonode)" | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;
# done


# # PFAM-Interpro connections
# for r in 1 2 3 4;
# do 
# 	echo "call apoc.load.json(\"file:///ontology/data/pfam-interpro/part$r.json\") yield value as entry
# 	with entry.metadata as datum
# 	with datum where datum.integrated is not null
# 	merge (inode:InterProFamily{family_id: datum.integrated})
# 	merge (pnode:PFAMFamily{family_id: datum.accession, family_type:datum.type})
# 	merge (inode)-[:interpro_pfam]-(pnode);" | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;

# done


# # Small Subunit nomenclature map
# echo "call apoc.load.json(\"file:///ontology/data/SSUMap.json\") yield value
# unwind(keys(value)) as key
# merge (pc:ProteinClass {class_id:key});" | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;

# # Large Subunit nomenclature map
# echo "\
# call apoc.load.json(\"file:///ontology/data/LSUMap.json\") yield value
# unwind(keys(value)) as key
# merge (pc:ProteinClass {class_id:key})" | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;

# # RNA nomenclature classes
# echo 'UNWIND [ 
#   "5SrRNA"  ,
#   "5.8SrRNA",
#   "12SrRNA" ,
#   "16SrRNA" ,
#   "21SrRNA" ,
#   "23SrRNA" ,
#   "25SrRNA" ,
#   "28SrRNA" ,
#   "35SrRNA" ,
#   "mRNA"    ,
#   "tRNA"    ]  as rnaclass
#   merge (n:RNAClass {class_id:rnaclass})' | cypher-shell --format plain --database $RIBOXYZ_DB_NAME;


