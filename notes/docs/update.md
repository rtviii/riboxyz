# The lifecycle of a single structure's ETL update

I consider a structure "processed" when all its mandatory assets (listed in `Assets` class) exist and its profile has been inducted into the database. 

## Assets:

- [M] profile 
- [M] mmcif strucuture
- PTC 
- exit tunnel
- thumbnail


## Assets Summary


<details>
  <summary>Assets Summary</summary>

### Profile 

The central part of the pipeline is the `ETLCollector` class which produces a `profile`. 
The gist of what it does:
    - query rcsb for various semantic parts of the structure (the node itself, the polymers, the ligands)
    - classify polymers
    - process structure metadata into riboxyz format (infer some semantic properties of interest)
    
### CIF Structure

MMCIF structure is downloaded locally to `RIBETL_DATA` and is used to generate structural landmarks etc.


### PTC

PTC Coordinates are rendered from the few reference structures


### Exit Tunnel is rendered  -- TODO

Presupposes PTC, render via mesh scripts

</details>
