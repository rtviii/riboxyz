## Background of the work:

- Define ribosome tunnel:
- Goals
- Khanh's 2019 NAR work (which led to riboxyz)
- Anton's evolutionary perspective

## Protocol:

### Point Cloud transformations:

- Mole method for extracting the tunnel (describe Mole briefly)
- Coordinate extraction <- Centerline [PIC]
- Coordinate extraction <- Bbox [PIC]
- Coordinate sphere expansion <- Bbox [PIC]
- Translate the whole thing to origin (index grid)
- Invert to get negative space

### Voxel Grid interior tunnel space capture:

Describe the problem of capturing the interior space of the tunnel, complications.
Describe DBSCAN briefly, its parameters, assumption of variable density.
Main cluster extraction:
- get cluster with dbscan
    (show that that the method is robust to parameter changes (show negative))
- Tunnel Point Cloud -> get surface points via delaunay triangulation
- Estimate and correct normals on the surface pointset
- Apply poisson reconstruction to the surface pointset


# TODOs:

- Clip via ribosome surface or the exit "crown" (get rid of the bulky solvent parts)
- Get PTC via RNA landmarks


