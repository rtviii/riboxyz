
# Extract `ribctl` from `api` structure

- [x]  there ought to be two separate clis for both

# Datastructure

- [x] separate assemblies
- [x] split ligand classes into instances 
- [ ] fit ligand parsing, transposition and fileformats to the new "polymeric factors" schema

# Structure Processing 

- [ ] integrate hmm-based classification into the process pipeline
- [ ] expand classification to rnas, trnas, factors
- [ ] improve ligand recognition ( BIRD, CCD?)


- [ ] rewrite taxonomy inference:
    - src ids and host ids separately
    - do toplevel struct even have a src_id in PDB?????
    - use "taxonomy.tax_infer_by_proportions" for a first-hand notion of what the structure is (add a new field ot the schema)

- [ ] flesh out the taxonomy : [src + host] [nodes in the db, hierarchy in the app(pull in the lib, construct dynamically on app start from assets)]


# Media

- [ ] render images for new structs 
- [ ] re-render profiles with the new ligand schema and separate assemblies 
- [ ] re-render ligands for all structs

- ligands/bsites:
    - render :
        - endpoint
        - typed [x]
        - logs
    - predict: 
        - endpoint
        - typed [x]
        - logs

- superimpose :
    - typed [x] 
    - endpoint 
    - logs

- [ ] RAII for peripheral data: PTC viewpoints, new nomenclature, tunnel obstruction profiles,

# Debugging 

- [ ] visualization troubleshoot

# DevOps:

- one-click docker deployment 
- how can we benefit from github actions?

# Later:

- "last updated" section
- automatic update per struct -- cron job every 24 weeks

- [x] rewrite struct init scripts in python
- [x] rewrtie cli in python

# LIGANDS

- ligands/ligandlike should have a common category in the db. filter the ions too. some ligands are not redndered at all (ex. 5AFI.FME)

- gene ontology + bird + CCD


## Static Files Server

Only after moving to UBC-ARC

- https://docs.djangoproject.com/en/4.1/howto/static-files/deployment/


## Integrate the PTC/conserved sites scripts scripts into the processing pipeline

