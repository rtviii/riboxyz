pull latest version

```sh
export NEO4J_URI=bolt://localhost:7687
export NEO4J_PASSWORD=... #usually empty
export NEO4J_USER=... # usually neo4j
export NEO4J_CURRENTDB=... # usually neo4j

export SECRET_KEY=... # django key
export DEBUG=False # django debug
export DJANGO_STATIC_FILES_PATH=... # django static files path, usually /api/staticfiles

# where the neo4j data is mounted
# can do /opt, but i just dump everything intot the user ease of use
export NEO4J_MOUNTPATH=/home/rtviii/dev/ribxz_neo4j_volume 

export RIBETL_DATA=... #main data path
export RIBXZ_ROOT=... # library code, usually'/home/$(whomai)/dev/riboxyz'
export CHIMERAX_BINARY_PATH=... # chimerax binary ex. '/usr/bin/chimerax'
```

# priamry assets/ribetl:

`p3 ribctl/ribd.py etl sync_all -t MMCIF -t STRUCTURE_PROFILE`
 
# neo4j docker setup:

if entirely new instance:
    
`p3 ribctl/ribd.py db instance create`
`p3 ribctl/ribd.py db upload-all`

