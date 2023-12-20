

# Consolidation tomorrow:

- tally each type of polymer, mito/cyto structures
- get rid of the garbage
- auto-update functionality
- boot up django, fix mixing docker-compose dependencies

<!-- ---------- -->
7p7q ->spectinomycin

spectonomycin -> 4wu1
spectonomycin -> 4wro
spectonomycin -> 6o8z


# WIP


### Datastructure:

- [ ] Pydantic 2.0

- [x] rewrite taxonomy inference: **<<<<<<<<<<**;
    - src ids and host ids separately
    - use "taxonomy.tax_infer_by_proportions" for a first-hand notion of what the structure is (add a new field ot the schema)
    - [x] flesh out the taxonomy : [src + host] [nodes in the db, hierarchy in the app(pull in the lib, construct dynamically on app start from assets)]

- [ ] expand classification to rnas, trnas, factors
    - [ ] integrate hmm-based classification into the process pipeline:
        - [ ] rRNA
        - [x] Proteins
        - [ ] Initiation Factors
        - [ ] Elongation Factors
        - [ ] Recycling Factors
        - [ ] tRNA


- [ ] Ligand improvements ( you can assume one has medicinal value/signficance if it's in drugbank)
    ### TODO:
    - [x] Added drugbank and pubchem references, that should be enough for now
    - [ ] Make sure the new types are propagated everywhere (including the api)
    - [ ] fit ligand parsing, transposition and fileformats to the new "polymeric factors" schema



# Roadmap

### Extract `ribctl` from `api` structure

- [x]  there ought to be two separate clis for both

### Datastructure

- [x] separate assemblies
- [x] split ligand classes into instances


### Media

- [ ] render images for new structs
- [ ] re-render profiles with the new ligand schema and separate assemblies
- [ ] re-render ligands for all structs

- ligands/bsites:
    - render:
        - endpoint
        - typed [x]
        - logs
    - predict:
        - endpoint
        - typed [x]
        - logs

- superimpose:
    - typed [x]
    - endpoint
    - logs

# Debugging

- [ ] visualization troubleshoot

# DevOps:

- one-click docker deployment
- how can we benefit from github actions?

# Later:

- "last updated" section
- automatic update per struct -- cron job every 24 weeks


# LIGANDS

- ligands/ligandlike should have a common category in the db. filter the ions too. some ligands are not redndered at all (ex. 5AFI.FME)

- gene ontology + bird + CCD



## Quirky parts

### Pymol

Pymol, so long as it serves its purpose for filling in the gaps of visual and geometric processing, should be installed the following way (taken verbatim from `docker-compose` file):

```bash
COPY pymol_source $PYMOL_SOURCE
ADD pymol_source $PYMOL_SOURCE
ENV PYMOL_PATH="${PYMOL_SOURCE}/__pymol_lib"
ENV PYTHONPATH="${PYMOL_PATH}/modules" 
RUN mkdir -p $PYMOL_PATH

WORKDIR $PYMOL_SOURCE
RUN python3 setup.py build install --home="${PYMOL_PATH}" --install-lib="${PYMOL_PATH}/modules/" --install-scripts="${PYMOL_PATH}"
# somewhere in the ~/.XXXrc's:
export PYTHONPATH="${PYTHONPATH}:${PYMOL_PATH}/modules"
```
Looking at the above by parts:
- OSS pymol installation requires the mmtf cpp packages. We grab those as is and they are already inside `pymol_source/include`.

- we install pymol headlessly throguh `setup.py` specifying a particular $PYMOL_PATH (a path of my choosing, it's only important that this ends up on pythonpath)

- to use pymol as a library, what we are really interested in is `$PYMOL_PATH/module` dir. We add that to $PYTHONPATH and voila, `from pymol import cmd` works anywhere.


The following pathing is elaborated on here: https://sourceforge.net/p/pymol/mailman/message/35988916/