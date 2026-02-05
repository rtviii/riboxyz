Alrighty, i have an NPET extraction pipeline in my riboxyz codebase that i want to clean up slightly (bubble up the main computational parameters).

 I already started by migrating the main parts of the codebase to "npet2" and it seems mostly done with the exception of some methods..

I have since implemented some caching, better run registration/manifest and refactored the whole thing into semi sensible stages. The major logic addition is also that i added a "grid refinement" step where basically we decrease the size of the voxel grid as second pass in the hopes of getting better separation and finer clusters/atom detail for occupied regions and hence more fidelity in the tunnel space.

My issue right now is the naming of artifacts and displaying of clusters in my visualization scrpt.. I also would like to produce a mesh from the bigges dbscan cluster at both grid sizes (1 and 0.5) to see how they compare and whether we gain something as well as some easier mechanism to play with eps/minnbrs params for each stage. 

These things i'd like to be able to visualize -- for each grid size: both the final mesh, the clusters of empty space and the bounding box. I mean this already works in my npet2_enhanced_viewer for the first stage, i think we just have to expand the namespcing a tad bit to accommodate different sets of clusters for each grid size and sets of eps/minnbrs parameters. I think we have a somehwat flexible manifest system going so i trust you to do this correctly. Also right now things get hashed in the "runs" directory into a fairly nondescript hash. let's at least append a date to it so the most recent bubbles up (we can truncate the date when checking whether a run with the given parameters exist already)...

Also there is a bunch of old prototyping parameters in the config.py right now -- if you see some that are definitely not getting used any where -- you can delete them so they dont' clutter up the picture.

Ok let me show you the code now. I should say that `lib/npet` is the old implemnetation (im happy to share more files from there, just ask me ) and `lib/npet2` is the new one that i'm refactoring out. 

Here is the codebase..


```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.npy|*.ply|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale' .  e      [npet_refactor]
.
├── __scripts
│   ├── gettrna.py
│   ├── ligclass_distill_composite.py
│   ├── ligclass_distill.py
│   ├── ligclass.py
│   ├── maps_acquire.py
│   ├── masif_plugin.py
│   ├── merge_lig_info.py
│   ├── mesh_grid.py
│   ├── move_classification.sh
│   ├── move_tunnels.sh
│   ├── narstructs_tally
│   ├── npet2_slice_viewer.py
│   ├── npet2_view.py
│   ├── npet2_viz_enhanced.py
│   ├── npet2_viz_run.py
│   ├── npy_to_pdb.py
│   ├── process_prediction_7k00.py
│   ├── pymol_visualtion.py
│   ├── rcsb_cumulative_entries_barplot.py
│   ├── reclassify.py
│   ├── ribxz_chimerax
│   │   ├── build
│   │   │   ├── bdist.macosx-10.9-universal2
│   │   │   └── lib
│   │   │       └── ribxz_chimerax
│   │   │           ├── __init__.py
│   │   │           ├── cmd_loci.py
│   │   │           ├── cmd_polymers.py
│   │   │           ├── cmd_registry.py
│   │   │           └── io.py
│   │   ├── bundle_info.xml
│   │   ├── dist
│   │   │   └── ribxz_chimerax-0.1-py3-none-any.whl
│   │   ├── ribxz_chimerax.egg-info
│   │   │   ├── dependency_links.txt
│   │   │   ├── PKG-INFO
│   │   │   ├── requires.txt
│   │   │   ├── SOURCES.txt
│   │   │   └── top_level.txt
│   │   └── src
│   │       ├── __init__.py
│   │       ├── cmd_loci.py
│   │       ├── cmd_polymers.py
│   │       ├── cmd_registry.py
│   │       └── io.py
│   ├── stats.json
│   ├── uniprot_seeds_query_record.py
│   └── williamson_assembly.py
├── dirtocontext.py
├── docker-compose.yml
├── Dockerfile-django
├── docs.md
├── muscle3.8.31_src.tar.gz
├── ncbi_taxonomy.sqlite
├── neo4j_ribosome
│   ├── __archive
│   │   ├── cypher_ops
│   │   │   ├── cypher_exec
│   │   │   ├── neo4j_commit_structure.sh
│   │   │   └── neo4j_seed_db_ontology.sh
│   │   └── riboxyz_seed_data
│   │       ├── ban-pfam-map-lsu.json
│   │       ├── ban-pfam-map-ssu.json
│   │       ├── interpro-base.json
│   │       ├── interpro-go-1.json
│   │       ├── interpro-go-2.json
│   │       ├── interpro-go-3.json
│   │       ├── interpro-go-4.json
│   │       ├── package-lock.json
│   │       ├── package.json
│   │       ├── pfam-interpro-1.json
│   │       ├── pfam-interpro-2.json
│   │       ├── pfam-interpro-3.json
│   │       ├── pfam-interpro-4.json
│   │       └── pfam-to-interpro-map.json
│   ├── __init__.py
│   ├── cypher
│   │   └── list_filtered_structs.cypher
│   ├── db_driver.py
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── node_ligand.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_protein.py
│   ├── node_rna.py
│   └── node_structure.py
├── neo4j.conf
├── neo4j.old.conf
├── notes
│   ├── binding_affinities.md
│   ├── docs
│   │   ├── architecture_environment_configs.md
│   │   └── update.md
│   ├── factors.md
│   ├── functional_sites.md
│   ├── general_directions.md
│   ├── hmm-based-classification.md
│   └── pymol.md
├── npet_orchestrator.py
├── NPET2
│   └── runs
│       └── 7K00
│           ├── 22cb6c9ebe9377ed
│           │   ├── manifest.json
│           │   └── stage
│           │       ├── 00_inputs
│           │       ├── 10_landmarks
│           │       ├── 20_exterior_shell
│           │       ├── 30_region_atoms
│           │       ├── 40_empty_space
│           │       ├── 50_clustering
│           │       ├── 55_grid_refine
│           │       ├── 60_surface_normals
│           │       └── 70_mesh_validate
│           ├── 27afaf798c882a2d
│           │   ├── manifest.json
│           │   └── stage
│           │       ├── 00_inputs
│           │       ├── 10_landmarks
│           │       ├── 20_exterior_shell
│           │       ├── 30_region_atoms
│           │       ├── 40_empty_space
│           │       ├── 50_clustering
│           │       ├── 55_grid_refine
│           │       ├── 60_surface_normals
│           │       └── 70_mesh_validate
│           ├── 759a0e4f9d385162
│           │   ├── manifest.json
│           │   └── stage
│           │       ├── 00_inputs
│           │       ├── 10_landmarks
│           │       ├── 20_exterior_shell
│           │       ├── 30_region_atoms
│           │       ├── 40_empty_space
│           │       ├── 50_clustering
│           │       ├── 60_surface_normals
│           │       └── 70_mesh_validate
│           └── 8c54c85531b076be
│               ├── manifest.json
│               └── stage
│                   ├── 00_inputs
│                   ├── 10_landmarks
│                   ├── 20_exterior_shell
│                   ├── 30_region_atoms
│                   ├── 40_empty_space
│                   ├── 50_clustering
│                   ├── 55_grid_refine
│                   ├── 60_surface_normals
│                   └── 70_mesh_validate
├── npet2_viewer_usage_examples.md
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_1.md
├── PLAN_refactor_npet_pipeline_2.md
├── PLAN_refactor_npet_pipeline_3.md
├── PLAN_refactor_npet_pipeline_4.md
├── PLAN_refactor_npet_pipeline_5.md
├── PLAN_refactor_npet_pipeline.md
├── q.md
├── ribctl
│   ├── __init__.py
│   ├── asset_manager
│   │   ├── asset_manager.py
│   │   ├── asset_registry.py
│   │   ├── asset_types.py
│   │   ├── doc.md
│   │   └── parallel_acquisition.py
│   ├── classifyre.json
│   ├── etl
│   │   ├── __init__.py
│   │   ├── etl_collector.py
│   │   └── gql_querystrings.py
│   ├── global_ops.py
│   ├── lib
│   │   ├── __libseq.py
│   │   ├── chimerax
│   │   │   ├── _cmd_ribrepr.py
│   │   │   ├── cmd_chainsplitter.py
│   │   │   ├── cmd_ligvis.py
│   │   │   ├── cmd_ribetl.py
│   │   │   ├── cmd_ribmovie.py
│   │   │   ├── cmd_ribrepr.py
│   │   │   ├── cmds_all.py
│   │   │   ├── ffmpeg_convert.sh
│   │   │   ├── ffmpeg_firstframe.sh
│   │   │   ├── gen_movies.py
│   │   │   ├── loop_ligvis.py
│   │   │   ├── loop_movies.py
│   │   │   ├── loop_split_chains.py
│   │   │   ├── notes.md
│   │   │   ├── produce_gif.py
│   │   │   └── thumbnails_from_ribetl.sh
│   │   ├── enumunion.py
│   │   ├── info.py
│   │   ├── landmarks
│   │   │   ├── __ptc_via_doris.py
│   │   │   ├── constriction_site.py
│   │   │   ├── notes.md
│   │   │   ├── ptc_via_trna.py
│   │   │   └── rrna_helices
│   │   │       ├── convert.py
│   │   │       ├── ecoli_7K00.json
│   │   │       ├── rrna_helices.py
│   │   │       ├── thermus_1VY4.json
│   │   │       └── yeast_4V88.json
│   │   ├── libbsite.py
│   │   ├── libhmm.py
│   │   ├── libmsa.py
│   │   ├── libseq.py
│   │   ├── libtax.py
│   │   ├── npet2
│   │   │   ├── __init__.py
│   │   │   ├── adapters
│   │   │   │   └── riboxyz_providers.py
│   │   │   ├── backends
│   │   │   │   ├── __init__.py
│   │   │   │   ├── clustering_io.py
│   │   │   │   ├── grid_occupancy.py
│   │   │   │   └── legacy
│   │   │   ├── core
│   │   │   │   ├── cache.py
│   │   │   │   ├── config.py
│   │   │   │   ├── interfaces.py
│   │   │   │   ├── manifest.py
│   │   │   │   ├── pipeline.py
│   │   │   │   ├── run_id.py
│   │   │   │   ├── settings.py
│   │   │   │   ├── store.py
│   │   │   │   └── types.py
│   │   │   ├── run.py
│   │   │   └── stages
│   │   │       ├── bootstrap.py
│   │   │       ├── grid_refine.py
│   │   │       └── legacy_minimal.py
│   │   ├── nsearch_gemmi.py
│   │   ├── ribosome_types
│   │   ├── schema
│   │   │   ├── __init__.py
│   │   │   ├── primitives.py
│   │   │   ├── types_binding_site.py
│   │   │   └── types_ribosome.py
│   │   ├── seq_project_many_to_one.py
│   │   ├── thumbnail.py
│   │   ├── types
│   │   │   └── polymer
│   │   │       ├── __init__.py
│   │   │       ├── base.py
│   │   │       ├── hierarchies.py
│   │   │       └── types.py
│   │   └── utils.py
│   ├── logger_config.py
│   ├── logs
│   │   ├── etl.log
│   │   └── loggers.py
│   ├── ribd.py
│   └── ribosome_ops.py
├── taxdump.tar.gz
├── test_npet2.py
├── v_9GRqXU
├── v_HqJuBQ
├── v_lKRFXT
└── v_WKwp5u
e  [error opening dir]

81 directories, 183 files

```

ribctl/lib/npet2/adapters/riboxyz_providers.py
```py
# ribctl/lib/npet2/adapters/riboxyz_providers.py
from __future__ import annotations

from typing import Any, Dict
import numpy as np

from ribctl.ribosome_ops import RibosomeOps
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.landmarks.constriction_site import get_constriction


class RiboxyzStructureProvider:
    def fingerprint(self, rcsb_id: str) -> str:
        # You can improve later: checksum mmcif, assembly ID, etc.
        p = AssetType.MMCIF.get_path(rcsb_id)
        return f"mmcif:{p}"

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        ro = RibosomeOps(rcsb_id)
        structure = ro.assets.biopython_structure()
        # simplest: extract all atom coords for first model
        atoms = [a for a in structure[0].get_atoms()]
        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray([getattr(a, "element", "") or a.get_id()[0] for a in atoms])
        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(AssetType.MMCIF.get_path(rcsb_id)),
            "profile": ro.profile,
            "ro": ro,  # keep around for legacy stages; core doesn’t require it
        }


class RiboxyzLandmarkProvider:
    def fingerprint(self, rcsb_id: str) -> str:
        # encode algorithm choices here later
        return "ptc_via_trna+constriction_site:v1"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        ptc = np.array(PTC_location(rcsb_id).location, dtype=np.float32)
        constr = np.array(get_constriction(rcsb_id), dtype=np.float32)
        return {"ptc_xyz": ptc, "constriction_xyz": constr}

```

ribctl/lib/npet2/backends/__init__.py
```py

```

ribctl/lib/npet2/backends/clustering_io.py
```py
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

def clusters_from_labels(points: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray]:
    clusters: dict[int, list[int]] = {}
    for i, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(i)
    out = {}
    for lab, idxs in clusters.items():
        out[lab] = points[np.asarray(idxs, dtype=np.int32)]
    return out

def write_dbscan_pass(out_dir: Path, *, prefix: str, points: np.ndarray, labels: np.ndarray) -> dict:
    out_dir = out_dir / prefix
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(out_dir / "points.npy", points.astype(np.float32))
    np.save(out_dir / "labels.npy", labels.astype(np.int32))

    clusters = clusters_from_labels(points, labels)
    index = {"prefix": prefix, "n_points": int(points.shape[0]), "clusters": []}

    for cid, pts in sorted(clusters.items(), key=lambda kv: kv[0]):
        p = out_dir / f"cluster_id{cid}.npy"
        np.save(p, pts.astype(np.float32))
        index["clusters"].append({"id": int(cid), "n": int(pts.shape[0]), "path": p.name})

    (out_dir / "index.json").write_text(json.dumps(index, indent=2))
    return index

```

ribctl/lib/npet2/backends/grid_occupancy.py
```py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class GridSpec:
    # origin and voxel size define world coords: x = origin + i*voxel
    origin: np.ndarray          # (3,)
    voxel_size: float
    shape: Tuple[int, int, int] # (nx, ny, nz)


def make_cylinder_grid(radius_A: float, height_A: float, voxel_A: float) -> GridSpec:
    """
    Canonical cylinder in C0:
      x in [-R, R], y in [-R, R], z in [0, H]
    """
    nx = int(np.floor((2 * radius_A) / voxel_A)) + 1
    ny = int(np.floor((2 * radius_A) / voxel_A)) + 1
    nz = int(np.floor(height_A / voxel_A)) + 1
    origin = np.array([-radius_A, -radius_A, 0.0], dtype=np.float32)
    return GridSpec(origin=origin, voxel_size=float(voxel_A), shape=(nx, ny, nz))


def grid_world_coords(grid: GridSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns axis coordinate arrays (x, y, z) for voxel centers along each axis.
    """
    ox, oy, oz = grid.origin
    nx, ny, nz = grid.shape
    v = grid.voxel_size
    x = ox + np.arange(nx, dtype=np.float32) * v
    y = oy + np.arange(ny, dtype=np.float32) * v
    z = oz + np.arange(nz, dtype=np.float32) * v
    return x, y, z


def cylinder_mask(grid: GridSpec, radius_A: float) -> np.ndarray:
    """
    Boolean mask of voxels inside cylinder radius (in C0).
    """
    x, y, _ = grid_world_coords(grid)
    X, Y = np.meshgrid(x, y, indexing="ij")
    inside = (X * X + Y * Y) <= (radius_A * radius_A)
    # broadcast across z
    return inside[:, :, None]


def points_to_occupied_seeds(points_c0: np.ndarray, grid: GridSpec) -> np.ndarray:
    """
    Convert points (C0 coords) to a sparse occupied grid of seeds at nearest voxels.
    """
    pts = np.asarray(points_c0, dtype=np.float32)
    v = grid.voxel_size
    origin = grid.origin

    ijk = np.floor((pts - origin[None, :]) / v + 0.5).astype(np.int32)
    nx, ny, nz = grid.shape

    valid = (
        (ijk[:, 0] >= 0) & (ijk[:, 0] < nx) &
        (ijk[:, 1] >= 0) & (ijk[:, 1] < ny) &
        (ijk[:, 2] >= 0) & (ijk[:, 2] < nz)
    )
    ijk = ijk[valid]
    occ = np.zeros(grid.shape, dtype=np.bool_)
    if ijk.shape[0] > 0:
        occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return occ


def occupancy_via_edt(points_c0: np.ndarray, grid: GridSpec, atom_radius_A: float) -> np.ndarray:
    """
    Occupancy grid: voxel is occupied if within atom_radius_A of any atom center.

    Steps:
      - seed occupied at nearest voxels
      - edt on ~occupied gives distance (in voxels) to nearest seed
      - threshold <= r_vox
    """
    seeds = points_to_occupied_seeds(points_c0, grid)

    # Distance (in voxels) from each voxel to nearest True in 'seeds'
    # distance_transform_edt computes distance to nearest zero;
    # so we compute on ~seeds, where zeros correspond to seeds.
    dist_vox = ndimage.distance_transform_edt(~seeds)

    r_vox = float(atom_radius_A) / float(grid.voxel_size)
    occupied = dist_vox <= r_vox
    return occupied


def empty_points_from_mask(grid: GridSpec, empty_mask: np.ndarray) -> np.ndarray:
    """
    Return coordinates (C0) of voxel centers where empty_mask is True.
    """
    empty_idx = np.where(empty_mask)
    x, y, z = grid_world_coords(grid)
    pts = np.column_stack((x[empty_idx[0]], y[empty_idx[1]], z[empty_idx[2]])).astype(np.float32)
    return pts


def save_grid_npy(grid: GridSpec, data: np.ndarray, path: Path, *, compress: bool = False) -> None:
    """
    Save a 3D grid along with its GridSpec metadata.
    File format: {path}_data.npy + {path}_spec.json
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if compress:
        np.savez_compressed(data_path.with_suffix('.npz'), data=data)
    else:
        np.save(data_path, data)
    
    spec_dict = {
        "origin": grid.origin.tolist(),
        "voxel_size": float(grid.voxel_size),
        "shape": list(grid.shape),
    }
    spec_path.write_text(__import__('json').dumps(spec_dict, indent=2))


def load_grid_npy(path: Path) -> tuple[GridSpec, np.ndarray]:
    """
    Load a grid saved by save_grid_npy.
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if data_path.with_suffix('.npz').exists():
        data = np.load(data_path.with_suffix('.npz'))['data']
    else:
        data = np.load(data_path)
    
    spec_dict = __import__('json').loads(spec_path.read_text())
    grid = GridSpec(
        origin=np.array(spec_dict["origin"], dtype=np.float32),
        voxel_size=float(spec_dict["voxel_size"]),
        shape=tuple(spec_dict["shape"]),
    )
    return grid, data


def voxel_to_world(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    """
    Convert voxel indices (i,j,k) to world coordinates.
    ijk: (N, 3) or (3,) array of voxel indices
    Returns: (N, 3) or (3,) world coordinates
    """
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin + ijk * grid.voxel_size


def world_to_voxel(grid: GridSpec, xyz: np.ndarray) -> np.ndarray:
    """
    Convert world coordinates to voxel indices.
    xyz: (N, 3) or (3,) world coordinates
    Returns: (N, 3) or (3,) voxel indices (floats; use floor/round as needed)
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    return (xyz - grid.origin) / grid.voxel_size


def get_occupied_voxel_centers(grid: GridSpec, occupancy: np.ndarray) -> np.ndarray:
    """
    Get world coordinates of occupied voxel centers.
    """
    occupied_idx = np.argwhere(occupancy)
    return voxel_to_world(grid, occupied_idx)
# ribctl/lib/npet2/backends/grid_occupancy.py
# ADD these functions at the end:

from scipy import ndimage

def connected_components_3d(
    binary_mask: np.ndarray, 
    connectivity: int = 26
) -> tuple[np.ndarray, int]:
    """
    Find connected components in a 3D binary mask.
    
    Args:
        binary_mask: 3D boolean array
        connectivity: 6 (face), 18 (face+edge), or 26 (face+edge+corner)
    
    Returns:
        labeled: Array same shape as input with component labels (0=background)
        n_components: Number of components found
    """
    if connectivity == 6:
        structure = ndimage.generate_binary_structure(3, 1)
    elif connectivity == 18:
        structure = ndimage.generate_binary_structure(3, 2)
    elif connectivity == 26:
        structure = ndimage.generate_binary_structure(3, 3)
    else:
        raise ValueError(f"connectivity must be 6, 18, or 26, got {connectivity}")
    
    labeled, n_components = ndimage.label(binary_mask, structure=structure)
    return labeled, n_components


def get_largest_component(labeled: np.ndarray, n_components: int) -> np.ndarray:
    """
    Extract mask of the largest connected component.
    
    Args:
        labeled: Output from connected_components_3d
        n_components: Number of components
    
    Returns:
        Binary mask of largest component only
    """
    if n_components == 0:
        return np.zeros_like(labeled, dtype=bool)
    
    # Count voxels in each component (excluding background=0)
    component_sizes = np.bincount(labeled.ravel())
    component_sizes[0] = 0  # Ignore background
    
    largest_label = np.argmax(component_sizes)
    return labeled == largest_label


def get_component_stats(labeled: np.ndarray, n_components: int) -> list[dict]:
    """
    Get statistics for all connected components.
    
    Returns:
        List of dicts with {label, size, bbox_min, bbox_max}
    """
    stats = []
    
    for label in range(1, n_components + 1):
        mask = labeled == label
        size = int(mask.sum())
        
        if size == 0:
            continue
        
        indices = np.argwhere(mask)
        bbox_min = indices.min(axis=0)
        bbox_max = indices.max(axis=0)
        
        stats.append({
            "label": int(label),
            "size": size,
            "bbox_min": bbox_min.tolist(),
            "bbox_max": bbox_max.tolist(),
        })
    
    # Sort by size descending
    stats.sort(key=lambda x: x["size"], reverse=True)
    return stats


def morphological_clean(
    binary_mask: np.ndarray, 
    operation: str = "opening",
    iterations: int = 1
) -> np.ndarray:
    """
    Apply morphological operations to clean up a binary mask.
    
    Args:
        binary_mask: 3D boolean array
        operation: "opening" (remove small bits), "closing" (fill small holes), 
                  "erosion", "dilation"
        iterations: Number of times to apply operation
    
    Returns:
        Cleaned binary mask
    """
    if operation == "opening":
        return ndimage.binary_opening(binary_mask, iterations=iterations)
    elif operation == "closing":
        return ndimage.binary_closing(binary_mask, iterations=iterations)
    elif operation == "erosion":
        return ndimage.binary_erosion(binary_mask, iterations=iterations)
    elif operation == "dilation":
        return ndimage.binary_dilation(binary_mask, iterations=iterations)
    else:
        raise ValueError(f"Unknown operation: {operation}")
```

ribctl/lib/npet2/core/cache.py
```py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import shutil

from .run_id import stable_hash_dict

@dataclass(frozen=True)
class StageCacheKey:
    stage: str
    inputs_fp: dict
    params: dict
    impl_version: str = "v1"  # bump when you change semantics

    def digest(self) -> str:
        return stable_hash_dict({
            "stage": self.stage,
            "inputs_fp": self.inputs_fp,
            "params": self.params,
            "impl_version": self.impl_version,
        })[:20]

class LocalStageCache:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def entry_dir(self, key: StageCacheKey) -> Path:
        d = self.root / key.stage / key.digest()
        d.mkdir(parents=True, exist_ok=True)
        return d

    def has(self, key: StageCacheKey, required: list[str]) -> bool:
        d = self.root / key.stage / key.digest()
        if not d.exists():
            return False
        return all((d / r).exists() for r in required)

    def copy_into(self, key: StageCacheKey, dest: Path, files: list[str]) -> None:
        src = self.root / key.stage / key.digest()
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            shutil.copy2(src / f, dest / f)

    def put_from(self, key: StageCacheKey, src_dir: Path, files: list[str]) -> None:
        dst = self.entry_dir(key)
        for f in files:
            shutil.copy2(src_dir / f, dst / f)

```

ribctl/lib/npet2/core/config.py
```py
# ribctl/lib/npet2/core/config.py
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Literal, Optional


@dataclass(frozen=True)
class DBSCANPolicy:
    eps_A: float                     # physical eps in Å
    density_factor: float = 0.10     # α in your back-of-envelope
    min_samples_override: Optional[int] = None


from dataclasses import dataclass, field
from typing import List, Literal


@dataclass(frozen=True)
class GridLevelConfig:
    name: str
    voxel_size_A: float

    atom_radius_mode: Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0

    occupancy_backend: Literal["legacy_kdtree", "grid_stamp", "edt", "gpu"] = "legacy_kdtree"
    roi_backend: Literal["full_cylinder", "bbox_from_prev", "tube_from_prev"] = "full_cylinder"


@dataclass(frozen=True)
class RunConfig:
    # region definition
    cylinder_radius_A: float = 35.0
    cylinder_height_A: float = 120.0

    # exterior shell (legacy alpha stage parameters)
    alpha_d3d_alpha: float = 200
    alpha_d3d_tol: float = 10
    alpha_d3d_offset: float = 3
    alpha_kdtree_radius: float = 40
    alpha_max_nn: int = 60
    alpha_tangent_planes_k: int = 20
    alpha_poisson_depth: int = 6
    alpha_poisson_ptweight: int = 4
    alpha_fill_holes: float = 2000

    # clustering (legacy values)
    dbscan_eps_A: float = 5.5
    dbscan_min_samples: int = 600
    refine_eps_A: float = 3.5
    refine_min_samples: int = 175

    # surface extraction
    surface_alpha: float = 2
    surface_tolerance: float = 1
    surface_offset: float = 2

    # normals
    normals_radius: float = 10
    normals_max_nn: int = 15
    normals_tangent_k: int = 10

    # mesh reconstruction
    mesh_poisson_depth: int = 6
    mesh_poisson_ptweight: int = 3

    # --- grid refinement (0.5Å ROI pass)
    refine_voxel_size_A      : float = 0.5
    refine_roi_pad_A         : float = 10.0
    refine_atom_radius_A     : float = 2.0
    refine_cc_connectivity   : int   = 26
    refine_topk_preview      : int   = 5
    refine_max_preview_points: int   = 50_000

    # grid refinement robustness
    refine_occ_close_iters: int = 1        # seals 0.5Å cracks
    refine_keep_within_A: float = 8.0      # keep void near coarse tunnel (0 disables)
    refine_forbid_roi_boundary: bool = True
    # optionally bump atom radius a bit for fine pass:
    refine_atom_radius_A: float = 2.5

    refine_voxel_size_A = 0.5
    refine_keep_within_A = 6.0      # keep local, but not too tight (try 5–8)
    refine_void_open_iters = 1      # helps break bridges
    refine_occ_close_iters = 0

    # DBSCAN on boundary points
    refine_dbscan_coarse_eps_A = 3
    refine_dbscan_coarse_min_samples = 30

    refine_dbscan_refine_eps_A = 3
    refine_dbscan_refine_min_samples = 20

    # if it’s too slow:
    # refine_dbscan_max_points = 250000
    # refine_dbscan_seed = 0

    # refinement plan (kept for later; default only one level)
    grid_levels: List[GridLevelConfig] = field(default_factory=lambda: [
        GridLevelConfig(name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"),
    ])


```

ribctl/lib/npet2/core/interfaces.py
```py
# ribctl/lib/npet2/core/interfaces.py
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Protocol, Tuple

import numpy as np

from .types import ArtifactRef, ArtifactType


class StructureProvider(Protocol):
    """
    Minimal structure access. Implemented by riboxyz adapters.
    """

    def fingerprint(self, rcsb_id: str) -> str:
        ...

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        """
        Return at minimum:
          - atom_xyz: (N,3) float32
          - atom_element: (N,) optional
        Can include:
          - mmcif_path, assemblies, chain ids, etc.
        """
        ...


class LandmarkProvider(Protocol):
    def fingerprint(self, rcsb_id: str) -> str:
        ...

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        """
        Must return:
          - ptc_xyz: (3,)
          - constriction_xyz: (3,)
        """
        ...


class ArtifactStore(Protocol):
    """
    Stores artifacts into the run directory and updates the manifest.
    """

    @property
    def run_dir(self) -> Path:
        ...

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType, data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        ...

    def put_json(self, *, name: str, stage: str, obj: Any, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        ...

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        ...

    def add_ref(self, ref: ArtifactRef) -> None:
        ...

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None:
        ...

```

ribctl/lib/npet2/core/manifest.py
```py
# ribctl/lib/npet2/core/manifest.py
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import time

from .types import ArtifactRef, ArtifactType


@dataclass
class StageRecord:
    name: str
    status: str = "pending"  # pending|running|success|failure|skipped
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    params: Dict[str, Any] = field(default_factory=dict)
    note: Optional[str] = None


@dataclass
class RunManifest:
    rcsb_id: str
    run_id: str
    pipeline_version: str
    created_at: float = field(default_factory=lambda: time.time())

    inputs: Dict[str, Any] = field(default_factory=dict)
    config_resolved: Dict[str, Any] = field(default_factory=dict)

    stages: Dict[str, StageRecord] = field(default_factory=dict)
    artifacts: List[Dict[str, Any]] = field(default_factory=list)

    success: Optional[bool] = None
    error: Optional[str] = None

    def add_artifact(self, ref: ArtifactRef) -> None:
        self.artifacts.append({
            "name": ref.name,
            "type": ref.type.value,
            "path": str(ref.path),
            "stage": ref.stage,
            "meta": ref.meta,
            "depends_on": list(ref.depends_on),
        })

    def to_json(self) -> str:
        # dataclasses → dict
        d = asdict(self)
        # StageRecord needs manual flatten
        d["stages"] = {k: asdict(v) for k, v in self.stages.items()}
        return json.dumps(d, indent=2)

    @staticmethod
    def from_path(path: Path) -> "RunManifest":
        data = json.loads(path.read_text())
        m = RunManifest(
            rcsb_id=data["rcsb_id"],
            run_id=data["run_id"],
            pipeline_version=data["pipeline_version"],
            created_at=data.get("created_at", time.time()),
            inputs=data.get("inputs", {}),
            config_resolved=data.get("config_resolved", {}),
        )
        m.success = data.get("success")
        m.error = data.get("error")
        # stages
        for k, v in data.get("stages", {}).items():
            m.stages[k] = StageRecord(**v)
        m.artifacts = data.get("artifacts", [])
        return m

```

ribctl/lib/npet2/core/pipeline.py
```py
# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List
import time

from .types import StageContext


class Stage(ABC):
    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]: ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None: ...


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        for stage in self.stages:
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)

            t0 = time.perf_counter()
            print(f"[npet2] >>> {stage.key} start")

            try:
                stage.run(ctx)
                dt = time.perf_counter() - t0
                print(f"[npet2] <<< {stage.key} done in {dt:,.2f}s")
                ctx.store.end_stage(stage.key, success=True, note=f"elapsed_s={dt:.3f}")
            except Exception as e:
                dt = time.perf_counter() - t0
                print(f"[npet2] !!! {stage.key} FAILED after {dt:,.2f}s: {e}")
                ctx.store.end_stage(
                    stage.key, success=False, note=f"elapsed_s={dt:.3f} err={e}"
                )
                ctx.store.finalize(success=False, error=str(e))
                raise

        ctx.store.finalize(success=True)
        return ctx

```

ribctl/lib/npet2/core/run_id.py
```py
# ribctl/lib/npet2/core/run_id.py
from __future__ import annotations

import hashlib
import json
from typing import Any, Dict


def stable_hash_dict(d: Dict[str, Any]) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def compute_run_id(*, rcsb_id: str, pipeline_version: str, inputs_fp: Dict[str, str], config_resolved: Dict[str, Any]) -> str:
    """
    run_id = sha256({rcsb_id, pipeline_version, inputs fingerprints, resolved config})
    """
    blob = {
        "rcsb_id"         : rcsb_id.upper(),
        "pipeline_version": pipeline_version,
        "inputs"          : dict(sorted(inputs_fp.items())),
        "config"          : config_resolved,
    }
    return stable_hash_dict(blob)[:16]  # short but sufficient; use full if you prefer

```

ribctl/lib/npet2/core/settings.py
```py
from pathlib import Path

NPET2_ROOT = Path("/Users/rtviii/dev/riboxyz/NPET2")
NPET2_RUNS_ROOT = NPET2_ROOT / "runs"

```

ribctl/lib/npet2/core/store.py
```py
# ribctl/lib/npet2/core/store.py
from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from .interfaces import ArtifactStore
from .manifest import RunManifest, StageRecord
from .types import ArtifactRef, ArtifactType


class LocalRunStore(ArtifactStore):
    def __init__(self, run_dir: Path, manifest: RunManifest):
        self._run_dir       = run_dir
        self._manifest      = manifest
        self._manifest_path = run_dir / "manifest.json"
        self._run_dir.mkdir(parents=True, exist_ok=True)
        self._write_manifest()

    def _abs(self, p: Path) -> Path:
        return p if p.is_absolute() else (self.run_dir / p)

    def _rel(self, p: Path) -> str:
        p = self._abs(p)
        try:
            return str(p.relative_to(self.run_dir))
        except ValueError:
            # Not under run_dir; fall back to absolute string (still tracked)
            return str(p)


    def register_file(
        self,
        *,
        name: str,
        stage: str,
        type: ArtifactType,
        path: Path,
        meta: Optional[Dict[str, Any]] = None,
        depends_on: tuple[str, ...] = (),
    ) -> ArtifactRef:
        ap = self._abs(path)

        ref = ArtifactRef(
            name=name,
            type=type,
            path=self._rel(ap),
            stage=stage,
            meta=meta or {},
            depends_on=depends_on,
        )
        self.add_ref(ref)
        return ref


    @property
    def run_dir(self) -> Path:
        return self._run_dir

    @property
    def manifest(self) -> RunManifest:
        return self._manifest

    def stage_dir(self, stage: str) -> Path:
        d = self._run_dir / "stage" / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def _write_manifest(self) -> None:
        self._manifest_path.write_text(self._manifest.to_json())

    def add_ref(self, ref: ArtifactRef) -> None:
        self._manifest.add_artifact(ref)
        self._write_manifest()

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType, data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        suffix = {
            ArtifactType.JSON: ".json",
            ArtifactType.NUMPY: ".npy",
            ArtifactType.PNG: ".png",
            ArtifactType.TXT: ".txt",
            ArtifactType.PLY_MESH: ".ply",
            ArtifactType.PLY_PCD: ".ply",
        }[type]
        out = self.stage_dir(stage) / f"{name}{suffix}"
        out.write_bytes(data)
        ref = ArtifactRef(name=name, type=type, path=out, stage=stage, meta=meta or {})
        self.add_ref(ref)
        return ref

    def put_json(self, *, name: str, stage: str, obj: Any, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        data = json.dumps(obj, indent=2).encode("utf-8")
        return self.put_bytes(name=name, stage=stage, type=ArtifactType.JSON, data=data, meta=meta)

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        out = self.stage_dir(stage) / f"{name}.npy"
        np.save(out, arr)
        ref = ArtifactRef(name=name, type=ArtifactType.NUMPY, path=out, stage=stage, meta=meta or {})
        self.add_ref(ref)
        return ref

    # Stage status helpers (optional but useful)
    def begin_stage(self, stage: str, params: Optional[Dict[str, Any]] = None) -> None:
        rec = self._manifest.stages.get(stage) or StageRecord(name=stage)
        rec.status = "running"
        rec.started_at = time.time()
        rec.params = params or {}
        self._manifest.stages[stage] = rec
        self._write_manifest()

    def end_stage(self, stage: str, success: bool, note: Optional[str] = None) -> None:
        rec = self._manifest.stages.get(stage) or StageRecord(name=stage)
        rec.status = "success" if success else "failure"
        rec.ended_at = time.time()
        rec.note = note
        self._manifest.stages[stage] = rec
        self._write_manifest()

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None:
        self._manifest.success = success
        self._manifest.error = error
        self._write_manifest()

```

ribctl/lib/npet2/core/types.py
```py
# ribctl/lib/npet2/core/types.py
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Mapping, Optional


class ArtifactType(str, Enum):
    JSON = "json"
    NUMPY = "npy"
    PLY_MESH = "ply_mesh"
    PLY_PCD = "ply_pcd"
    PNG = "png"
    TXT = "txt"


@dataclass(frozen=True)
class ArtifactRef:
    """
    A stable handle to an on-disk artifact, referenced from the manifest.
    """
    name: str                 # semantic name: "ptc", "empty_points_level_0"
    type: ArtifactType
    path: Path                # absolute or run-relative; store decides
    stage: str                # stage key: "10_landmarks"
    meta: Dict[str, Any] = field(default_factory=dict)
    depends_on: tuple[str, ...] = ()  # artifact names (or ids later)


@dataclass
class StageContext:
    """
    Shared context passed through the pipeline.
    - inputs: raw objects needed by compute (arrays, coords, providers, etc.)
    - artifacts: ArtifactRef registry for cross-stage access
    - stats: cheap summaries to help decisions (counts/bounds, etc.)
    """
    run_id: str
    rcsb_id: str
    config: Any  # RunConfig (kept Any to avoid import cycles)
    store: Any   # ArtifactStore

    inputs: Dict[str, Any] = field(default_factory=dict)
    artifacts: Dict[str, ArtifactRef] = field(default_factory=dict)
    stats: Dict[str, Any] = field(default_factory=dict)

    def require(self, key: str) -> Any:
        if key not in self.inputs:
            raise KeyError(f"Missing required input: {key}")
        return self.inputs[key]

    def require_artifact(self, name: str) -> ArtifactRef:
        if name not in self.artifacts:
            raise KeyError(f"Missing required artifact: {name}")
        return self.artifacts[name]

```

ribctl/lib/npet2/stages/bootstrap.py
```py
# ribctl/lib/npet2/stages/bootstrap.py
from __future__ import annotations

from dataclasses import asdict
from typing import Any, Dict

import numpy as np

from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext


class Stage00Inputs(Stage):
    key = "00_inputs"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        # include only config fields that actually affect this stage
        return {}

    def run(self, ctx: StageContext) -> None:
        structure_provider = ctx.require("structure_provider")
        data = structure_provider.load_atoms(ctx.rcsb_id)

        atom_xyz = np.asarray(data["atom_xyz"], dtype=np.float32)
        ctx.inputs["atom_xyz"] = atom_xyz
        ctx.inputs["atom_element"] = data.get("atom_element", None)

        # Keep adapter objects in ctx.inputs for now to support legacy backends later
        # (core won't *require* them; stages/backends can choose to use them)
        for k in ("mmcif_path", "profile", "ro"):
            if k in data:
                ctx.inputs[k] = data[k]

        # Save minimal artifact for debugging + provenance
        ctx.artifacts["atom_xyz"] = ctx.store.put_numpy(
            name="atom_xyz",
            stage=self.key,
            arr=atom_xyz,
            meta={"shape": list(atom_xyz.shape), "dtype": str(atom_xyz.dtype)},
        )

        # Useful stats to carry forward
        mins = atom_xyz.min(axis=0)
        maxs = atom_xyz.max(axis=0)
        ctx.stats["atom_bounds"] = {"min": mins.tolist(), "max": maxs.tolist()}
        ctx.stats["n_atoms"] = int(atom_xyz.shape[0])


class Stage10Landmarks(Stage):
    key = "10_landmarks"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        landmark_provider = ctx.require("landmark_provider")
        lm = landmark_provider.get_landmarks(ctx.rcsb_id)

        ptc = np.asarray(lm["ptc_xyz"], dtype=np.float32)
        constr = np.asarray(lm["constriction_xyz"], dtype=np.float32)

        ctx.inputs["ptc_xyz"] = ptc
        ctx.inputs["constriction_xyz"] = constr

        ctx.artifacts["ptc"] = ctx.store.put_json(
            name="ptc",
            stage=self.key,
            obj={"location": ptc.tolist()},
            meta={"units": "A"},
        )
        ctx.artifacts["constriction_site"] = ctx.store.put_json(
            name="constriction_site",
            stage=self.key,
            obj={"location": constr.tolist()},
            meta={"units": "A"},
        )

```

ribctl/lib/npet2/stages/grid_refine.py
```py
# ribctl/lib/npet2/stages/grid_refine.py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Tuple

import numpy as np
import pyvista as pv
from scipy import ndimage
import time


from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    transform_points_from_C0,
)

from ribctl.lib.npet2.backends.grid_occupancy import (
    GridSpec,
    occupancy_via_edt,
    connected_components_3d,
)


def _make_bbox_grid(lo: np.ndarray, hi: np.ndarray, voxel: float) -> GridSpec:
    lo = np.asarray(lo, dtype=np.float32)
    hi = np.asarray(hi, dtype=np.float32)
    voxel = float(voxel)

    span = hi - lo
    # +1 so both ends are representable
    shape = tuple((np.ceil(span / voxel).astype(np.int32) + 1).tolist())
    return GridSpec(origin=lo, voxel_size=voxel, shape=shape)


def _voxel_centers_from_indices(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin[None, :] + ijk * float(grid.voxel_size)


def _points_to_ijk(grid: GridSpec, pts_c0: np.ndarray) -> np.ndarray:
    """Nearest-voxel mapping for points in C0 -> ijk indices."""
    v = float(grid.voxel_size)
    ijk = np.floor((pts_c0 - grid.origin[None, :]) / v + 0.5).astype(np.int32)
    return ijk


def _valid_ijk(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    nx, ny, nz = grid.shape
    m = (
        (ijk[:, 0] >= 0)
        & (ijk[:, 0] < nx)
        & (ijk[:, 1] >= 0)
        & (ijk[:, 1] < ny)
        & (ijk[:, 2] >= 0)
        & (ijk[:, 2] < nz)
    )
    return m


def _topk_component_stats(labeled: np.ndarray, k: int = 6) -> list[dict]:
    sizes = np.bincount(labeled.ravel())
    if sizes.size == 0:
        return []
    sizes[0] = 0
    if sizes.sum() == 0:
        return []

    top_labels = np.argsort(sizes)[::-1][:k]
    out = []
    for lab in top_labels:
        if lab == 0 or sizes[lab] == 0:
            continue
        idx = np.argwhere(labeled == lab)
        bbox_min = idx.min(axis=0)
        bbox_max = idx.max(axis=0)
        out.append(
            {
                "label": int(lab),
                "size": int(sizes[lab]),
                "bbox_min": bbox_min.tolist(),
                "bbox_max": bbox_max.tolist(),
            }
        )
    return out

```

ribctl/lib/npet2/stages/legacy_minimal.py
```py
from __future__ import annotations
import json
from pathlib import Path
import time
from typing import Any, Dict, List, Tuple
import numpy as np
import pyvista as pv
import open3d as o3d

from ribctl.lib.npet2.backends.grid_occupancy import (
    connected_components_3d,
    occupancy_via_edt,
)
from ribctl.lib.npet2.core.cache import StageCacheKey
from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from scipy import ndimage

# Legacy helpers (keep pipeline operational)
from ribctl.lib.npet.alphalib import (
    cif_to_point_cloud,
    fast_normal_estimation,
    quick_surface_points,
    validate_mesh_pyvista,
)
from ribctl.lib.npet.kdtree_approach import (
    apply_poisson_reconstruction,
    ribosome_entities,
    filter_residues_parallel,
    transform_points_to_C0,
    transform_points_from_C0,
    create_point_cloud_mask,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    ptcloud_convex_hull_points,
    estimate_normals,
)
from ribctl.lib.npet2.stages.grid_refine import (
    _make_bbox_grid,
    _points_to_ijk,
    _topk_component_stats,
    _valid_ijk,
    _voxel_centers_from_indices,
)


def _tunnel_debris_chains(rcsb_id: str, ro, profile) -> List[str]:
    # your legacy hardcoded exclusions
    tunnel_debris = {
        "3J7Z": ["a", "7"],
        "5GAK": ["z"],
        "5NWY": ["s"],
        "7A5G": ["Y2"],
        "9F1D": ["BK"],
    }
    rcsb_id = rcsb_id.upper()
    skip = tunnel_debris.get(rcsb_id, []).copy()

    # mitochondrial mL45 (best-effort)
    if getattr(profile, "mitochondrial", False):
        try:
            chain = ro.get_poly_by_polyclass("mL45")
            if chain is not None:
                skip.append(chain.auth_asym_id)
        except Exception:
            pass
    return skip


class Stage20ExteriorShell(Stage):
    key = "20_exterior_shell"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "d3d_alpha": c.alpha_d3d_alpha,
            "d3d_tol": c.alpha_d3d_tol,
            "d3d_offset": c.alpha_d3d_offset,
            "kdtree_radius": c.alpha_kdtree_radius,
            "max_nn": c.alpha_max_nn,
            "tangent_k": c.alpha_tangent_planes_k,
            "poisson_depth": c.alpha_poisson_depth,
            "poisson_ptweight": c.alpha_poisson_ptweight,
            "fill_holes": c.alpha_fill_holes,
        }

    def run(self, ctx: StageContext) -> None:
        stage_cache = ctx.require("stage_cache")
        inputs_fp = ctx.require("inputs_fp")
        params = self.params(ctx)

        key = StageCacheKey(
            stage=self.key,
            inputs_fp={"structure": inputs_fp["structure"]},
            params=params,
            impl_version="v1",
        )

        stage_dir = ctx.store.stage_dir(self.key)
        cached_files = [
            "alpha_shell.ply",
            "alpha_shell_quality.json",
            "alpha_normals.ply",
            "alpha_surface_points.npy",
            "ribosome_ptcloud.npy",
        ]

        if stage_cache.has(
            key, required=["alpha_shell.ply", "alpha_shell_quality.json"]
        ):
            stage_cache.copy_into(key, stage_dir, cached_files)

            quality = json.loads((stage_dir / "alpha_shell_quality.json").read_text())
            ctx.inputs["alpha_shell_path"] = str(stage_dir / "alpha_shell.ply")
            ctx.inputs["alpha_shell_watertight"] = bool(
                quality.get("watertight", False)
            )

            # register artifacts (paths now exist under stage_dir)
            ctx.store.register_file(
                name="alpha_shell_mesh",
                stage=self.key,
                type=ArtifactType.PLY_MESH,
                path=stage_dir / "alpha_shell.ply",
            )
            ctx.store.register_file(
                name="alpha_shell_quality",
                stage=self.key,
                type=ArtifactType.JSON,
                path=stage_dir / "alpha_shell_quality.json",
            )
            return

        c = ctx.config
        ro = ctx.require("ro")
        cifpath = Path(ctx.require("mmcif_path"))

        stage_dir = ctx.store.stage_dir(self.key)
        ptcloud_path = stage_dir / "ribosome_ptcloud.npy"
        surface_pts_path = stage_dir / "alpha_surface_points.npy"
        normals_pcd_path = stage_dir / "alpha_normals.ply"
        mesh_path = stage_dir / "alpha_shell.ply"
        quality_path = stage_dir / "alpha_shell_quality.json"

        # point cloud from cif (legacy)
        first_assembly_chains = ro.first_assembly_auth_asym_ids()
        ptcloud = cif_to_point_cloud(
            str(cifpath), first_assembly_chains, do_atoms=True
        ).astype(np.float32)
        np.save(ptcloud_path, ptcloud)
        ctx.store.register_file(
            name="ribosome_ptcloud",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=ptcloud_path,
        )

        # surface points
        surface_pts = quick_surface_points(
            ptcloud, c.alpha_d3d_alpha, c.alpha_d3d_tol, c.alpha_d3d_offset
        ).astype(np.float32)
        np.save(surface_pts_path, surface_pts)
        ctx.store.register_file(
            name="alpha_surface_points",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=surface_pts_path,
        )

        # normal estimation (legacy)
        normal_estimated_pcd = fast_normal_estimation(
            surface_pts, c.alpha_kdtree_radius, c.alpha_max_nn, c.alpha_tangent_planes_k
        )

        # robust-ish normal orientation: outward
        center = normal_estimated_pcd.get_center()
        normal_estimated_pcd.orient_normals_towards_camera_location(
            camera_location=center
        )
        normal_estimated_pcd.normals = o3d.utility.Vector3dVector(
            -np.asarray(normal_estimated_pcd.normals)
        )

        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        ctx.store.register_file(
            name="alpha_normals_pcd",
            stage=self.key,
            type=ArtifactType.PLY_PCD,
            path=normals_pcd_path,
        )

        # poisson reconstruction (writes mesh_path)
        apply_poisson_reconstruction(
            str(normals_pcd_path),
            mesh_path,
            recon_depth=c.alpha_poisson_depth,
            recon_pt_weight=c.alpha_poisson_ptweight,
        )

        # repair + keep largest component
        mesh = pv.read(mesh_path)
        mesh = mesh.fill_holes(c.alpha_fill_holes)
        mesh = mesh.connectivity(largest=True).triangulate()
        mesh.save(mesh_path)

        watertight = validate_mesh_pyvista(mesh)

        # record quality
        quality = {
            "watertight": bool(watertight),
            "n_points": int(mesh.n_points),
            "n_faces": int(mesh.n_faces),
            "open_edges": int(mesh.n_open_edges),
            "is_manifold": bool(mesh.is_manifold),
            "bounds": list(mesh.bounds),
        }
        quality_path.write_text(__import__("json").dumps(quality, indent=2))
        ctx.store.register_file(
            name="alpha_shell_quality",
            stage=self.key,
            type=ArtifactType.JSON,
            path=quality_path,
        )

        ctx.store.register_file(
            name="alpha_shell_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
        )
        ctx.inputs["alpha_shell_path"] = str(mesh_path)
        ctx.inputs["alpha_shell_watertight"] = bool(watertight)
        if watertight:
            stage_cache.put_from(key, stage_dir, cached_files)


class Stage30RegionAtoms(Stage):
    key = "30_region_atoms"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {"radius_A": c.cylinder_radius_A, "height_A": c.cylinder_height_A}

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        ro = ctx.require("ro")
        profile = ctx.require("profile")
        cifpath = ctx.require("mmcif_path")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        skip = _tunnel_debris_chains(ctx.rcsb_id, ro, profile)

        residues = ribosome_entities(
            rcsb_id=ctx.rcsb_id,
            cifpath=cifpath,
            level="R",
            skip_nascent_chain=skip,
        )

        filtered_residues = filter_residues_parallel(
            residues=residues,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,  # <- stops process fan-out
            chunk_size=5000,  # <- optional; irrelevant if max_workers=1
        )

        filtered_points = np.asarray(
            [atom.get_coord() for r in filtered_residues for atom in r.child_list],
            dtype=np.float32,
        )

        out = ctx.store.stage_dir(self.key) / "region_atom_xyz.npy"
        np.save(out, filtered_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out,
            meta={"n": int(filtered_points.shape[0])},
        )

        ctx.inputs["region_atom_xyz"] = filtered_points


class Stage40EmptySpace(Stage):
    key = "40_empty_space"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "grid_levels": [
                {
                    "name": gl.name,
                    "voxel_size_A": gl.voxel_size_A,
                    "backend": gl.occupancy_backend,
                    "atom_radius_mode": getattr(gl, "atom_radius_mode", "uniform"),
                    "uniform_atom_radius_A": getattr(gl, "uniform_atom_radius_A", None),
                }
                for gl in c.grid_levels
            ],
            "radius_A": c.cylinder_radius_A,
            "height_A": c.cylinder_height_A,
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import numpy as np
        import pyvista as pv

        from ribctl.lib.npet.kdtree_approach import (
            create_point_cloud_mask,
            transform_points_from_C0,
            transform_points_to_C0,
        )

        # EDT / grid backend
        from ribctl.lib.npet2.backends.grid_occupancy import (
            make_cylinder_grid,
            cylinder_mask,
            occupancy_via_edt,
            empty_points_from_mask,
            save_grid_npy,
            get_occupied_voxel_centers,
        )

        c = ctx.config
        region_xyz = np.asarray(ctx.require("region_atom_xyz"), dtype=np.float32)
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        # Transform once (invariant across levels)
        region_c0 = transform_points_to_C0(region_xyz, ptc, constr)

        # Load shell once
        shell = pv.read(alpha_shell_path)
        if not isinstance(shell, pv.PolyData):
            shell = shell.extract_surface()
        if not shell.is_all_triangles:
            shell = shell.triangulate()

        # If Stage20 recorded watertightness, use it to pick a safer mode
        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        stage_dir = ctx.store.stage_dir(self.key)

        # Optional: record a tiny stage note artifact about clipping mode
        clip_note = {
            "alpha_shell_path": str(alpha_shell_path),
            "alpha_shell_watertight": watertight,
            "clipping_mode": "select_enclosed_points(check_surface=True)"
            if watertight
            else "select_enclosed_points(check_surface=False) [fallback; shell not watertight]",
        }
        p_clip_note = stage_dir / "clipping_note.json"
        p_clip_note.write_text(json.dumps(clip_note, indent=2))
        ctx.store.register_file(
            name="clipping_note",
            stage=self.key,
            type=ArtifactType.JSON,
            path=p_clip_note,
        )

        last_empty = None

        for gl in c.grid_levels:
            backend = gl.occupancy_backend

            if backend == "legacy_kdtree":
                # Legacy semantics:
                # create_point_cloud_mask returns final_mask where:
                #   True = occupied OR outside-cylinder
                # so ~mask are empty voxels INSIDE the cylinder.
                mask, (x, y, z) = create_point_cloud_mask(
                    region_c0,
                    radius=c.cylinder_radius_A,
                    height=c.cylinder_height_A,
                    voxel_size=gl.voxel_size_A,
                    radius_around_point=gl.uniform_atom_radius_A,
                )

                idx = np.where(~mask)
                empty_c0 = np.column_stack((x[idx[0]], y[idx[1]], z[idx[2]])).astype(
                    np.float32
                )

            elif backend == "edt":
                # Grid / EDT backend:
                # - build canonical cylinder grid in C0
                # - compute occupancy within atom radius via EDT
                # - IMPORTANT: outside-cylinder must be treated as occupied (matches legacy)
                grid = make_cylinder_grid(
                    radius_A=float(c.cylinder_radius_A),
                    height_A=float(c.cylinder_height_A),
                    voxel_A=float(gl.voxel_size_A),
                )

                occ = occupancy_via_edt(
                    region_c0,
                    grid,
                    atom_radius_A=float(gl.uniform_atom_radius_A),
                )

                cyl2d = cylinder_mask(
                    grid, radius_A=float(c.cylinder_radius_A)
                )  # (nx,ny,1)
                cyl = np.broadcast_to(cyl2d, grid.shape)  # (nx,ny,nz)

                # Match legacy: outside cylinder is "occupied" so it never becomes empty.
                occ = occ | (~cyl)

                empty_mask = ~occ
                empty_c0 = empty_points_from_mask(grid, empty_mask & cyl)

                # Persist grids for visualization/debug
                # This writes:
                #   occupancy_grid_<name>_data.npy + occupancy_grid_<name>_spec.json
                #   empty_mask_<name>_data.npy     + empty_mask_<name>_spec.json
                save_grid_npy(
                    grid, occ, stage_dir / f"occupancy_grid_{gl.name}", compress=False
                )
                save_grid_npy(
                    grid,
                    (empty_mask & cyl),
                    stage_dir / f"empty_mask_{gl.name}",
                    compress=False,
                )

                # Register the grid files explicitly
                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"occupancy_grid_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"occupancy_grid_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"empty_mask_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"empty_mask_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )

                # Optional convenience artifact: occupied voxel centers (WORLD) for fast viz
                # (Your viz script already has logic for occupied_voxels_level_0.npy.)
                try:
                    occ_centers_c0 = get_occupied_voxel_centers(grid, occ).astype(
                        np.float32
                    )
                    occ_centers_world = transform_points_from_C0(
                        occ_centers_c0, ptc, constr
                    ).astype(np.float32)
                    p_occ = stage_dir / f"occupied_voxels_{gl.name}.npy"
                    np.save(p_occ, occ_centers_world)
                    ctx.store.register_file(
                        name=f"occupied_voxels_{gl.name}",
                        stage=self.key,
                        type=ArtifactType.NUMPY,
                        path=p_occ,
                        meta={
                            "voxel_size_A": gl.voxel_size_A,
                            "n": int(occ_centers_world.shape[0]),
                        },
                    )
                except Exception:
                    # keep EDT path robust even if convenience write fails
                    pass

            else:
                raise ValueError(
                    f"Grid level {gl.name}: unsupported backend {backend} "
                    f"(supported: legacy_kdtree, edt)"
                )

            # Convert empty voxel centers back to world coords
            empty_world = transform_points_from_C0(empty_c0, ptc, constr).astype(
                np.float32
            )

            # Optional: persist pre-clip points (helps debug non-watertight shells)
            p_pre = stage_dir / f"empty_points_{gl.name}_preclip.npy"
            np.save(p_pre, empty_world)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}_preclip",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=p_pre,
                meta={"voxel_size_A": gl.voxel_size_A, "n": int(empty_world.shape[0])},
            )

            # Clip to alpha shell interior (best effort)
            # If shell isn't watertight, check_surface=False avoids hard failure but may be imperfect.
            if empty_world.shape[0] == 0:
                inside = empty_world
            else:
                pts_poly = pv.PolyData(empty_world)
                sel = pts_poly.select_enclosed_points(shell, check_surface=watertight)
                inside = empty_world[sel["SelectedPoints"] == 1].astype(np.float32)

            # Save artifacts per level
            out = stage_dir / f"empty_points_{gl.name}.npy"
            np.save(out, inside)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "voxel_size_A": gl.voxel_size_A,
                    "n": int(inside.shape[0]),
                    "backend": backend,
                    "alpha_shell_watertight": watertight,
                },
            )

            ctx.inputs[f"empty_points_{gl.name}"] = inside
            last_empty = inside

        ctx.inputs["empty_points"] = last_empty


class Stage50Clustering(Stage):
    key = "50_clustering"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "dbscan_eps_A": c.dbscan_eps_A,
            "dbscan_min_samples": c.dbscan_min_samples,
            "refine_eps_A": c.refine_eps_A,
            "refine_min_samples": c.refine_min_samples,
            # toggles (safe defaults; you can bubble to config later)
            "save_all_clusters": bool(getattr(c, "dbscan_save_all_clusters", True)),
            "save_noise_cluster": bool(getattr(c, "dbscan_save_noise_cluster", False)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import time

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        # ---- inputs
        if "empty_points" not in ctx.inputs:
            # helpful debug for you if Stage40 didn’t set it
            print(f"[50_clustering] ctx.inputs keys: {sorted(ctx.inputs.keys())}")
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)
        if empty_pts.ndim != 2 or empty_pts.shape[1] != 3 or empty_pts.shape[0] == 0:
            raise ValueError(
                f"[50_clustering] empty_points must be (N,3) and non-empty, got {empty_pts.shape}"
            )

        save_all = bool(getattr(c, "dbscan_save_all_clusters", True))
        save_noise = bool(getattr(c, "dbscan_save_noise_cluster", False))

        print(f"[50_clustering] empty_points n={empty_pts.shape[0]:,}")

        # Helper: save a DBSCAN pass into subdir (points.npy, labels.npy, plus per-cluster files optionally)
        def _save_pass(
            pass_name: str,
            pts: np.ndarray,
            labels: np.ndarray,
            clusters_dict: dict[int, list],
            eps: float,
            min_samples: int,
        ) -> None:
            pdir = stage_dir / pass_name
            pdir.mkdir(parents=True, exist_ok=True)

            np.save(pdir / "points.npy", pts.astype(np.float32))
            np.save(pdir / "labels.npy", labels.astype(np.int32))

            # summary/index for quick browsing + viz logic
            counts = {}
            for lab in np.unique(labels):
                counts[int(lab)] = int((labels == lab).sum())
            index = {
                "pass": pass_name,
                "eps_A": float(eps),
                "min_samples": int(min_samples),
                "n_points": int(pts.shape[0]),
                "labels": counts,  # includes -1
            }
            (pdir / "index.json").write_text(json.dumps(index, indent=2))

            if not save_all:
                return

            # save per-cluster arrays
            # clusters_dict is label -> list[point]
            for lab, plist in clusters_dict.items():
                lab = int(lab)
                if lab == -1 and not save_noise:
                    continue
                arr = np.asarray(plist, dtype=np.float32)
                if arr.size == 0:
                    continue
                np.save(pdir / f"cluster_id{lab}.npy", arr)

        # -----------------------
        # PASS 1: coarse DBSCAN
        # -----------------------
        t0 = time.perf_counter()
        db0, clusters0 = DBSCAN_capture(empty_pts, c.dbscan_eps_A, c.dbscan_min_samples)
        labels0 = np.asarray(db0.labels_, dtype=np.int32)
        dt0 = time.perf_counter() - t0
        print(
            f"[50_clustering] coarse DBSCAN took {dt0:,.2f}s labels={len(set(labels0.tolist())):,}"
        )

        _save_pass(
            "coarse",
            empty_pts,
            labels0,
            clusters0,
            eps=c.dbscan_eps_A,
            min_samples=c.dbscan_min_samples,
        )

        # pick winner from coarse
        largest, largest_id = DBSCAN_pick_largest_cluster(clusters0)
        largest = np.asarray(largest, dtype=np.float32)
        if largest.ndim != 2 or largest.shape[1] != 3 or largest.shape[0] == 0:
            raise ValueError("[50_clustering] largest cluster is empty/unexpected")

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest)
        ctx.store.register_file(
            name="largest_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_largest,
            meta={"cluster_id": int(largest_id), "n": int(largest.shape[0])},
        )

        # -----------------------
        # PASS 2: refine DBSCAN
        # -----------------------
        t1 = time.perf_counter()
        db1, clusters1 = DBSCAN_capture(largest, c.refine_eps_A, c.refine_min_samples)
        labels1 = np.asarray(db1.labels_, dtype=np.int32)
        dt1 = time.perf_counter() - t1
        print(
            f"[50_clustering] refine DBSCAN took {dt1:,.2f}s labels={len(set(labels1.tolist())):,}"
        )

        _save_pass(
            "refine",
            largest,
            labels1,
            clusters1,
            eps=c.refine_eps_A,
            min_samples=c.refine_min_samples,
        )

        refined, refined_id = DBSCAN_pick_largest_cluster(clusters1)
        refined = np.asarray(refined, dtype=np.float32)
        if refined.ndim != 2 or refined.shape[1] != 3 or refined.shape[0] == 0:
            raise ValueError("[50_clustering] refined cluster is empty/unexpected")

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0])},
        )

        # ---- outputs for downstream
        ctx.inputs["refined_cluster"] = refined
        # optional: keep also the coarse winner handy
        ctx.inputs["largest_cluster"] = largest

        print(
            f"[50_clustering] winner coarse={largest.shape[0]:,} refine={refined.shape[0]:,}"
        )


class Stage55GridRefine05(Stage):
    """
    Stage55 (refined grid) using DBSCAN segmentation (like Stage50), not connected-components selection.

    Outputs (for viz, Stage60+):
      stage/55_grid_refine/
        refined_surface_points_level_1.npy      (WORLD coords; used downstream as refined_cluster)
        dbscan_coarse/points.npy + labels.npy   (WORLD coords; labels from DBSCAN in C0)
        dbscan_refine/points.npy + labels.npy   (WORLD coords; labels from DBSCAN in C0)
        dbscan_diagnostics.json
        roi_bbox_c0.json
        grid_spec_level_1.json
        void_mask_level_1.npy   (uint8, optional debug)
        occupied_mask_level_1.npy (uint8, optional debug)
    """

    key = "55_grid_refine"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A": float(getattr(c, "refine_voxel_size_A", 0.5)),
            "roi_pad_A": float(getattr(c, "refine_roi_pad_A", 10.0)),
            "atom_radius_A": float(getattr(c, "refine_atom_radius_A", 2.0)),

            # optional localization: restrict void to within this distance of Stage50 refined cluster (0=off)
            "keep_within_A": float(getattr(c, "refine_keep_within_A", 0.0)),

            # optional morphology (helps break skinny bridges BEFORE DBSCAN)
            "occ_close_iters": int(getattr(c, "refine_occ_close_iters", 0)),   # closing on occupied
            "void_open_iters": int(getattr(c, "refine_void_open_iters", 0)),   # opening on void (break bridges)

            # prevent planar ROI faces from polluting void
            "forbid_roi_boundary": bool(getattr(c, "refine_forbid_roi_boundary", True)),

            # DBSCAN params (coarse then refine)
            "dbscan_coarse_eps_A": float(getattr(c, "refine_dbscan_coarse_eps_A", 1.5)),
            "dbscan_coarse_min_samples": int(getattr(c, "refine_dbscan_coarse_min_samples", 30)),
            "dbscan_refine_eps_A": float(getattr(c, "refine_dbscan_refine_eps_A", 1.0)),
            "dbscan_refine_min_samples": int(getattr(c, "refine_dbscan_refine_min_samples", 20)),

            # safety caps (0 = no cap)
            "dbscan_max_points": int(getattr(c, "refine_dbscan_max_points", 0)),
            "dbscan_seed": int(getattr(c, "refine_dbscan_seed", 0)),

            # diagnostics
            "max_cluster_stats": int(getattr(c, "refine_dbscan_max_cluster_stats", 25)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        from pathlib import Path

        import numpy as np
        import pyvista as pv
        from scipy import ndimage
        from scipy.spatial import cKDTree
        from sklearn.cluster import DBSCAN

        c = ctx.config
        stage_dir = Path(ctx.store.stage_dir(self.key))
        stage_dir.mkdir(parents=True, exist_ok=True)

        # ----------------------------
        # inputs
        # ----------------------------
        refined_world = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)
        if refined_world.ndim != 2 or refined_world.shape[1] != 3 or refined_world.shape[0] == 0:
            raise ValueError(f"[55_grid_refine] refined_cluster invalid: {refined_world.shape}")

        region_xyz = np.asarray(ctx.require("region_atom_xyz"), dtype=np.float32)
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = str(ctx.require("alpha_shell_path"))
        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        voxel = float(getattr(c, "refine_voxel_size_A", 0.5))
        pad = float(getattr(c, "refine_roi_pad_A", 10.0))
        atom_r = float(getattr(c, "refine_atom_radius_A", 2.0))
        keep_within_A = float(getattr(c, "refine_keep_within_A", 0.0))
        occ_close_iters = int(getattr(c, "refine_occ_close_iters", 0))
        void_open_iters = int(getattr(c, "refine_void_open_iters", 0))
        forbid_roi_boundary = bool(getattr(c, "refine_forbid_roi_boundary", True))

        eps_c = float(getattr(c, "refine_dbscan_coarse_eps_A", 1.5))
        ms_c = int(getattr(c, "refine_dbscan_coarse_min_samples", 30))
        eps_r = float(getattr(c, "refine_dbscan_refine_eps_A", 1.0))
        ms_r = int(getattr(c, "refine_dbscan_refine_min_samples", 20))

        dbscan_max_points = int(getattr(c, "refine_dbscan_max_points", 0))
        dbscan_seed = int(getattr(c, "refine_dbscan_seed", 0))
        max_stats = int(getattr(c, "refine_dbscan_max_cluster_stats", 25))

        # ----------------------------
        # ROI in C0
        # ----------------------------
        refined_c0 = transform_points_to_C0(refined_world, ptc, constr).astype(np.float32)
        lo = refined_c0.min(axis=0)
        hi = refined_c0.max(axis=0)
        lo_pad = lo - pad
        hi_pad = hi + pad

        roi_obj = {
            "roi_id": "bbox_pad_dbscan",
            "frame": "C0",
            "pad_A": float(pad),
            "lo": [float(x) for x in lo_pad.tolist()],
            "hi": [float(x) for x in hi_pad.tolist()],
            "transform": {"ptc": [float(x) for x in ptc.tolist()], "constriction": [float(x) for x in constr.tolist()]},
            "source": {"stage": "50_clustering", "artifact": "refined_cluster"},
        }
        (stage_dir / "roi_bbox_c0.json").write_text(json.dumps(roi_obj, indent=2))
        ctx.artifacts["roi_bbox_c0"] = str(stage_dir / "roi_bbox_c0.json")

        # ----------------------------
        # select atoms near ROI in C0
        # ----------------------------
        region_c0 = transform_points_to_C0(region_xyz, ptc, constr).astype(np.float32)
        lo_sel = lo_pad - atom_r
        hi_sel = hi_pad + atom_r
        m_atoms = np.all((region_c0 >= lo_sel[None, :]) & (region_c0 <= hi_sel[None, :]), axis=1)
        atoms_roi_c0 = region_c0[m_atoms]
        if atoms_roi_c0.shape[0] == 0:
            raise ValueError("[55_grid_refine] no atoms selected near ROI")

        # ----------------------------
        # build bbox grid + cylinder mask
        # ----------------------------
        grid = _make_bbox_grid(lo_pad, hi_pad, voxel)

        def _cylinder_mask_bbox_grid(grid: GridSpec, radius_A: float, zmin_A: float, zmax_A: float) -> np.ndarray:
            nx, ny, nz = grid.shape
            v = float(grid.voxel_size)
            ox, oy, oz = grid.origin

            x = ox + np.arange(nx, dtype=np.float32) * v
            y = oy + np.arange(ny, dtype=np.float32) * v
            z = oz + np.arange(nz, dtype=np.float32) * v

            X, Y = np.meshgrid(x, y, indexing="ij")
            inside_r = (X * X + Y * Y) <= (radius_A * radius_A)
            inside_z = (z >= zmin_A) & (z <= zmax_A)
            return inside_r[:, :, None] & inside_z[None, None, :]

        cyl = _cylinder_mask_bbox_grid(
            grid,
            radius_A=float(c.cylinder_radius_A),
            zmin_A=0.0,
            zmax_A=float(c.cylinder_height_A),
        )

        # ----------------------------
        # occupancy via EDT
        # ----------------------------
        occupied = occupancy_via_edt(atoms_roi_c0, grid, atom_radius_A=atom_r)

        if occ_close_iters > 0:
            occupied = ndimage.binary_closing(occupied, iterations=occ_close_iters)

        # outside cylinder is "occupied" so void can't leak
        occupied = occupied | (~cyl)

        np.save(stage_dir / "occupied_mask_level_1.npy", occupied.astype(np.uint8))

        # ----------------------------
        # empty points (C0) inside cylinder
        # ----------------------------
        empty_mask = (~occupied) & cyl
        empty_idx = np.argwhere(empty_mask)
        if empty_idx.shape[0] == 0:
            raise ValueError("[55_grid_refine] empty_mask had 0 voxels in ROI")

        empty_pts_c0 = _voxel_centers_from_indices(grid, empty_idx).astype(np.float32)

        # ----------------------------
        # clip empties to inside alpha shell interior (C0)
        # ----------------------------
        shell_world = pv.read(alpha_shell_path)
        if not isinstance(shell_world, pv.PolyData):
            shell_world = shell_world.extract_surface()
        shell_world = shell_world.triangulate()

        shell_c0 = shell_world.copy(deep=True)
        shell_c0.points = transform_points_to_C0(np.asarray(shell_world.points, dtype=np.float32), ptc, constr)

        sel = pv.PolyData(empty_pts_c0).select_enclosed_points(shell_c0, check_surface=watertight)
        inside_flags = (np.asarray(sel["SelectedPoints"], dtype=np.int8) == 1)

        inside_idx = empty_idx[inside_flags]
        void_mask = np.zeros_like(empty_mask, dtype=np.bool_)
        if inside_idx.shape[0] == 0:
            raise ValueError("[55_grid_refine] no empty voxels inside alpha shell in ROI")
        void_mask[inside_idx[:, 0], inside_idx[:, 1], inside_idx[:, 2]] = True

        # forbid ROI boundary to avoid planar faces contaminating the void
        if forbid_roi_boundary:
            void_mask[0, :, :] = False
            void_mask[-1, :, :] = False
            void_mask[:, 0, :] = False
            void_mask[:, -1, :] = False
            void_mask[:, :, 0] = False
            void_mask[:, :, -1] = False

        # optional localization around Stage50 refined cluster
        if keep_within_A > 0.0:
            coarse_ijk = _points_to_ijk(grid, refined_c0)
            m_valid = _valid_ijk(grid, coarse_ijk)
            coarse_ijk = coarse_ijk[m_valid]
            if coarse_ijk.shape[0] > 0:
                seed = np.zeros(grid.shape, dtype=np.bool_)
                seed[coarse_ijk[:, 0], coarse_ijk[:, 1], coarse_ijk[:, 2]] = True
                dist_vox = ndimage.distance_transform_edt(~seed)
                r_vox = float(keep_within_A) / float(voxel)
                void_mask = void_mask & (dist_vox <= r_vox)

        # optional opening on void to break skinny bridges
        if void_open_iters > 0:
            st = ndimage.generate_binary_structure(3, 1)  # 6-neighborhood
            void_mask = ndimage.binary_opening(void_mask, structure=st, iterations=void_open_iters)

        np.save(stage_dir / "void_mask_level_1.npy", void_mask.astype(np.uint8))

        # ----------------------------
        # boundary extraction (surface points) in C0
        # ----------------------------
        st_er = ndimage.generate_binary_structure(3, 1)
        er = ndimage.binary_erosion(void_mask, structure=st_er, iterations=1)
        boundary_mask = void_mask & (~er)
        boundary_idx = np.argwhere(boundary_mask)
        if boundary_idx.shape[0] == 0:
            raise ValueError("[55_grid_refine] boundary extraction produced 0 voxels")

        boundary_pts_c0 = _voxel_centers_from_indices(grid, boundary_idx).astype(np.float32)

        # ----------------------------
        # DBSCAN helper (score clusters by closeness to Stage50 refined tunnel in C0)
        # ----------------------------
        tree = cKDTree(refined_c0)

        def _cluster_stats(points_c0: np.ndarray, labels: np.ndarray) -> list[dict]:
            stats = []
            for lab in np.unique(labels):
                if lab == -1:
                    continue
                m = labels == lab
                pts = points_c0[m]
                if pts.shape[0] == 0:
                    continue
                d, _ = tree.query(pts, k=1)
                stats.append(
                    {
                        "label": int(lab),
                        "size": int(pts.shape[0]),
                        "median_dist_to_stage50_A": float(np.median(d)),
                        "p05_dist_to_stage50_A": float(np.percentile(d, 5)),
                        "p95_dist_to_stage50_A": float(np.percentile(d, 95)),
                    }
                )
            stats.sort(key=lambda x: (x["median_dist_to_stage50_A"], -x["size"]))
            return stats

        def _choose_best_label(stats: list[dict]) -> int:
            if not stats:
                return -1
            # prefer smallest median distance; tie-break by size
            return int(stats[0]["label"])

        def _maybe_cap(points_c0: np.ndarray, points_w: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
            """
            Optionally cap points for DBSCAN speed by random subsample.
            Returns: (points_c0_cap, points_w_cap, idx_cap)
            """
            n = points_c0.shape[0]
            if dbscan_max_points and n > dbscan_max_points:
                rng = np.random.default_rng(dbscan_seed)
                k = int(dbscan_max_points)
                idx = rng.choice(n, size=k, replace=False)
                idx.sort()
                return points_c0[idx], points_w[idx], idx
            idx = np.arange(n, dtype=np.int64)
            return points_c0, points_w, idx

        # boundary points in WORLD for saving/viz
        boundary_pts_w = transform_points_from_C0(boundary_pts_c0, ptc, constr).astype(np.float32)

        # cap if requested
        boundary_pts_c0_cap, boundary_pts_w_cap, idx_cap = _maybe_cap(boundary_pts_c0, boundary_pts_w)

        # ----------------------------
        # coarse DBSCAN on refined-grid boundary points
        # ----------------------------
        db_coarse = DBSCAN(eps=eps_c, min_samples=ms_c, metric="euclidean", n_jobs=-1)
        labels_coarse = db_coarse.fit_predict(boundary_pts_c0_cap).astype(np.int32)

        d_coarse = stage_dir / "dbscan_coarse"
        d_coarse.mkdir(parents=True, exist_ok=True)
        np.save(d_coarse / "points.npy", boundary_pts_w_cap.astype(np.float32))
        np.save(d_coarse / "labels.npy", labels_coarse.astype(np.int32))

        stats_coarse = _cluster_stats(boundary_pts_c0_cap, labels_coarse)
        best_coarse = _choose_best_label(stats_coarse)
        if best_coarse == -1:
            raise ValueError(
                f"[55_grid_refine] coarse DBSCAN produced no clusters (all noise). "
                f"Try increasing eps or decreasing min_samples. eps={eps_c} min_samples={ms_c}"
            )

        # points for refine pass
        m_best = labels_coarse == best_coarse
        pts_refine_c0 = boundary_pts_c0_cap[m_best]
        pts_refine_w = boundary_pts_w_cap[m_best]
        if pts_refine_c0.shape[0] == 0:
            raise ValueError("[55_grid_refine] best coarse cluster had 0 points (unexpected)")

        # ----------------------------
        # refine DBSCAN inside best coarse cluster
        # ----------------------------
        db_refine = DBSCAN(eps=eps_r, min_samples=ms_r, metric="euclidean", n_jobs=-1)
        labels_refine = db_refine.fit_predict(pts_refine_c0).astype(np.int32)

        d_ref = stage_dir / "dbscan_refine"
        d_ref.mkdir(parents=True, exist_ok=True)
        np.save(d_ref / "points.npy", pts_refine_w.astype(np.float32))
        np.save(d_ref / "labels.npy", labels_refine.astype(np.int32))

        stats_refine = _cluster_stats(pts_refine_c0, labels_refine)
        best_refine = _choose_best_label(stats_refine)
        if best_refine == -1:
            raise ValueError(
                f"[55_grid_refine] refine DBSCAN produced no clusters (all noise). "
                f"Try increasing eps or decreasing min_samples. eps={eps_r} min_samples={ms_r}"
            )

        final_mask = labels_refine == best_refine
        final_surface_w = pts_refine_w[final_mask].astype(np.float32)
        if final_surface_w.shape[0] == 0:
            raise ValueError("[55_grid_refine] final selected refine cluster had 0 points")

        # ----------------------------
        # diagnostics + grid spec
        # ----------------------------
        spec_obj = {
            "frame": "C0",
            "origin": [float(x) for x in grid.origin.tolist()],
            "voxel_size_A": float(grid.voxel_size),
            "shape": [int(x) for x in grid.shape],
            "transform": {"ptc": [float(x) for x in ptc.tolist()], "constriction": [float(x) for x in constr.tolist()]},
        }
        (stage_dir / "grid_spec_level_1.json").write_text(json.dumps(spec_obj, indent=2))
        ctx.inputs["grid_spec_level_1_path"] = str(stage_dir / "grid_spec_level_1.json")

        diag = {
            "voxel_size_A": voxel,
            "roi_pad_A": pad,
            "atom_radius_A": atom_r,
            "keep_within_A": keep_within_A,
            "occ_close_iters": occ_close_iters,
            "void_open_iters": void_open_iters,
            "forbid_roi_boundary": forbid_roi_boundary,
            "alpha_shell_watertight": watertight,
            "boundary_points_total": int(boundary_pts_c0.shape[0]),
            "boundary_points_used_for_dbscan": int(boundary_pts_c0_cap.shape[0]),
            "dbscan_max_points": int(dbscan_max_points),
            "dbscan_coarse": {
                "eps_A": eps_c,
                "min_samples": ms_c,
                "best_label": int(best_coarse),
                "n_clusters": int(len(stats_coarse)),
                "clusters": stats_coarse[:max_stats],
            },
            "dbscan_refine": {
                "eps_A": eps_r,
                "min_samples": ms_r,
                "best_label": int(best_refine),
                "n_clusters": int(len(stats_refine)),
                "clusters": stats_refine[:max_stats],
            },
            "final_surface_points": int(final_surface_w.shape[0]),
        }
        (stage_dir / "dbscan_diagnostics.json").write_text(json.dumps(diag, indent=2))

        print(
            f"[55_grid_refine] DBSCAN refined-grid boundary: "
            f"voxel={voxel} void_vox={int(void_mask.sum()):,} boundary_vox={int(boundary_mask.sum()):,} "
            f"coarse eps={eps_c} ms={ms_c} best={best_coarse} -> refine eps={eps_r} ms={ms_r} best={best_refine} "
            f"final_pts={final_surface_w.shape[0]:,} cap={dbscan_max_points}"
        )

        # ----------------------------
        # outputs for Stage60+
        # ----------------------------
        np.save(stage_dir / "refined_surface_points_level_1.npy", final_surface_w)

        ctx.inputs["refined_cluster_surface"] = True
        ctx.inputs["refined_cluster"] = final_surface_w

        # keep artifacts consistent with your viz script expectations
        ctx.artifacts["refined_surface_points_level_1"] = str(stage_dir / "refined_surface_points_level_1.npy")
        ctx.artifacts["refined_cluster_level_1"] = str(stage_dir / "refined_surface_points_level_1.npy")







class Stage60SurfaceNormals(Stage):
    key = "60_surface_normals"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "surface_alpha": c.surface_alpha,
            "surface_tolerance": c.surface_tolerance,
            "surface_offset": c.surface_offset,
            "normals_radius": c.normals_radius,
            "normals_max_nn": c.normals_max_nn,
            "normals_tangent_k": c.normals_tangent_k,
        }

    def run(self, ctx: StageContext) -> None:
        import time

        c = ctx.config
        refined = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)

        surface_flag = bool(ctx.inputs.get("refined_cluster_surface", False))
        print(
            f"[60_surface_normals] refined_cluster n={refined.shape[0]:,} surface_flag={surface_flag}"
        )

        stage_dir = ctx.store.stage_dir(self.key)

        # Decide how to get surface points
        if surface_flag:
            # Points already represent a surface-like sampling (e.g., boundary voxels at 0.5Å)
            surface_pts = refined
            print(
                "[60_surface_normals] using refined points directly as surface_pts (skip Delaunay)"
            )
        else:
            # Old behavior (expensive on large point clouds)
            t0 = time.perf_counter()
            surface_pts = ptcloud_convex_hull_points(
                refined, c.surface_alpha, c.surface_tolerance, c.surface_offset
            ).astype(np.float32)
            dt = time.perf_counter() - t0
            print(
                f"[60_surface_normals] delaunay_3d+extract_surface took {dt:,.2f}s surface_pts n={surface_pts.shape[0]:,}"
            )

        # Save surface points
        p_surface = stage_dir / "surface_points.npy"
        np.save(p_surface, surface_pts)
        ctx.store.register_file(
            name="surface_points",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_surface,
            meta={"n": int(surface_pts.shape[0])},
        )

        # Normals estimation
        t1 = time.perf_counter()
        pcd = estimate_normals(
            surface_pts,
            kdtree_radius=c.normals_radius,
            kdtree_max_nn=c.normals_max_nn,
            correction_tangent_planes_n=c.normals_tangent_k,
        )
        dt1 = time.perf_counter() - t1
        print(f"[60_surface_normals] estimate_normals took {dt1:,.2f}s")

        # Write normals point cloud
        p_normals = stage_dir / "surface_normals.ply"
        o3d.io.write_point_cloud(str(p_normals), pcd)
        ctx.store.register_file(
            name="surface_normals_pcd",
            stage=self.key,
            type=ArtifactType.PLY_PCD,
            path=p_normals,
        )

        ctx.inputs["normals_pcd_path"] = str(p_normals)


class Stage70MeshValidate(Stage):
    key = "70_mesh_validate"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "poisson_depth": c.mesh_poisson_depth,
            "poisson_ptweight": c.mesh_poisson_ptweight,
            # optional knobs for voxel-mesh cleanup
            "voxel_fill_holes_A": float(getattr(c, "voxel_mesh_fill_holes_A", 50.0)),
            "voxel_smooth_iters": int(getattr(c, "voxel_mesh_smooth_iters", 0)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import numpy as np
        import pyvista as pv

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        mesh_path = stage_dir / "npet2_tunnel_mesh.ply"

        # ---- helper: validate mesh with logging
        def _mesh_stats(m: pv.PolyData) -> dict:
            return {
                "n_points": int(m.n_points),
                "n_faces": int(m.n_faces),
                "open_edges": int(m.n_open_edges),
                "is_manifold": bool(m.is_manifold),
                "bounds": [float(x) for x in m.bounds],
            }

        # ---- 1) Try legacy Poisson route first (if we have normals)
        normals_pcd_path = ctx.inputs.get("normals_pcd_path", None)
        if normals_pcd_path:
            try:
                apply_poisson_reconstruction(
                    str(normals_pcd_path),
                    mesh_path,
                    recon_depth=c.mesh_poisson_depth,
                    recon_pt_weight=c.mesh_poisson_ptweight,
                )
            except Exception as e:
                print(f"[70_mesh_validate] poisson threw exception: {e}")

        # If Poisson produced a mesh, validate it
        if mesh_path.exists():
            try:
                m = pv.read(str(mesh_path))
                st = _mesh_stats(m)
                print(f"[70_mesh_validate] poisson mesh stats: {st}")
                watertight = validate_mesh_pyvista(m)
                if watertight:
                    ctx.store.register_file(
                        name="tunnel_mesh",
                        stage=self.key,
                        type=ArtifactType.PLY_MESH,
                        path=mesh_path,
                        meta={"watertight": True, "method": "poisson"},
                    )
                    ctx.inputs["tunnel_mesh_path"] = str(mesh_path)
                    return
                else:
                    print(
                        "[70_mesh_validate] poisson mesh not watertight; falling back to voxel meshing"
                    )
            except Exception as e:
                print(
                    f"[70_mesh_validate] failed reading/validating poisson mesh; falling back: {e}"
                )
        else:
            print(
                "[70_mesh_validate] poisson did not produce a mesh file; falling back to voxel meshing"
            )

        # ---- 2) Fallback: deterministic voxel meshing (marching cubes / contour)
        mask_p = ctx.inputs.get("selected_void_component_mask_level_1_path", None)
        spec_p = ctx.inputs.get("grid_spec_level_1_path", None)

        if not (mask_p and spec_p and Path(mask_p).exists() and Path(spec_p).exists()):
            raise ValueError(
                "Final mesh is not watertight and voxel fallback inputs are missing "
                "(expected selected_void_component_mask_level_1_path + grid_spec_level_1_path)"
            )

        spec = json.loads(Path(spec_p).read_text())
        voxel = float(spec["voxel_size_A"])
        origin = np.asarray(spec["origin"], dtype=np.float32)

        ptc = np.asarray(spec["transform"]["ptc"], dtype=np.float32)
        constr = np.asarray(spec["transform"]["constriction"], dtype=np.float32)

        vol = np.load(mask_p).astype(np.uint8)
        if vol.ndim != 3:
            raise ValueError(
                f"[70_mesh_validate] voxel volume must be 3D, got {vol.shape}"
            )

        # Pad volume by 1 voxel so surfaces at the ROI boundary get "capped" (prevents open surfaces).
        vol_pad = np.pad(vol, 1, constant_values=0)
        origin_pad = origin - voxel  # because we added a 1-voxel pad on the low side

        # Build VTK ImageData: note cell_data length must match (nx*ny*nz)
        img = pv.ImageData(
            dimensions=(
                vol_pad.shape[0] + 1,
                vol_pad.shape[1] + 1,
                vol_pad.shape[2] + 1,
            ),
            spacing=(voxel, voxel, voxel),
            origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
        )
        img.cell_data["void"] = vol_pad.ravel(order="F")

        # Extract surface at 0.5 (binary volume)
        surf_c0 = img.contour(isosurfaces=[0.5], scalars="void").triangulate()
        if surf_c0.n_points == 0 or surf_c0.n_faces == 0:
            raise ValueError("[70_mesh_validate] voxel contour produced empty surface")

        # Keep largest connected surface, clean
        surf_c0 = surf_c0.clean(tolerance=0.0).connectivity(largest=True)

        # Optional cleanup: fill small holes (units are approx in Angstroms; VTK wants a size in "mesh units")
        fill_holes_A = float(getattr(c, "voxel_mesh_fill_holes_A", 50.0))
        try:
            surf_c0 = surf_c0.fill_holes(fill_holes_A)
        except Exception:
            pass

        # Optional smoothing (usually keep off unless needed)
        smooth_iters = int(getattr(c, "voxel_mesh_smooth_iters", 0))
        if smooth_iters > 0:
            try:
                surf_c0 = surf_c0.smooth(n_iter=smooth_iters)
            except Exception:
                pass

        # Transform mesh points C0 -> world
        pts_c0 = np.asarray(surf_c0.points, dtype=np.float32)
        pts_w = transform_points_from_C0(pts_c0, ptc, constr).astype(np.float32)
        surf_w = surf_c0.copy(deep=True)
        surf_w.points = pts_w

        # Compute normals for consistency (not required for watertightness, but nice to have)
        try:
            surf_w = surf_w.compute_normals(
                auto_orient_normals=True, consistent_normals=True
            )
        except Exception:
            pass

        surf_w.save(str(mesh_path))

        # Validate
        st2 = _mesh_stats(surf_w)
        print(f"[70_mesh_validate] voxel mesh stats: {st2}")

        watertight = validate_mesh_pyvista(surf_w)
        if not watertight:
            raise ValueError(
                "Final mesh is not watertight (voxel fallback also failed)"
            )

        ctx.store.register_file(
            name="tunnel_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"watertight": True, "method": "voxel_contour"},
        )
        ctx.inputs["tunnel_mesh_path"] = str(mesh_path)

```

ribctl/lib/npet2/run.py
```py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

from ribctl.lib.npet2.adapters.riboxyz_providers import (
    RiboxyzLandmarkProvider,
    RiboxyzStructureProvider,
)
from ribctl.lib.npet2.core.config import RunConfig
from ribctl.lib.npet2.core.manifest import RunManifest
from ribctl.lib.npet2.core.run_id import compute_run_id
from ribctl.lib.npet2.core.settings import NPET2_ROOT, NPET2_RUNS_ROOT
from ribctl.lib.npet2.core.store import LocalRunStore
from ribctl.lib.npet2.core.types import StageContext
from ribctl.lib.npet2.core.pipeline import Pipeline


from ribctl.lib.npet2.stages.bootstrap import Stage00Inputs, Stage10Landmarks
from ribctl.lib.npet2.stages.legacy_minimal import (
    Stage20ExteriorShell,
    Stage30RegionAtoms,
    Stage40EmptySpace,
    Stage50Clustering,
    Stage55GridRefine05,
    Stage60SurfaceNormals,
    Stage70MeshValidate,
)


def _pipeline_version() -> str:
    return "npet2-dev"


def run_npet2(
    rcsb_id: str,
    config: Optional[RunConfig] = None,
    *,
    structure_provider=None,
    landmark_provider=None,
) -> StageContext:
    rcsb_id = rcsb_id.upper()
    config = config or RunConfig()


    structure_provider = structure_provider or RiboxyzStructureProvider()
    landmark_provider = landmark_provider or RiboxyzLandmarkProvider()

    config_resolved = asdict(config)
    inputs_fp = {
        "structure": structure_provider.fingerprint(rcsb_id),
        "landmarks": landmark_provider.fingerprint(rcsb_id),
    }

    run_id = compute_run_id(
        rcsb_id=rcsb_id,
        pipeline_version=_pipeline_version(),
        inputs_fp=inputs_fp,
        config_resolved=config_resolved,
    )

    run_dir = NPET2_RUNS_ROOT / rcsb_id / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    manifest = RunManifest(
        rcsb_id=rcsb_id,
        run_id=run_id,
        pipeline_version=_pipeline_version(),
        inputs={"fingerprints": inputs_fp},
        config_resolved=config_resolved,
    )
    store = LocalRunStore(run_dir=run_dir, manifest=manifest)

    ctx = StageContext(
        run_id=run_id,
        rcsb_id=rcsb_id,
        config=config,
        store=store,
        inputs={
            "structure_provider": structure_provider,
            "landmark_provider": landmark_provider,
        },
    )

    # in run_npet2()
    from ribctl.lib.npet2.core.cache import LocalStageCache
    ctx.inputs["stage_cache"] = LocalStageCache(NPET2_ROOT / "cache")
    ctx.inputs["inputs_fp"] = inputs_fp

    pipeline = Pipeline(
        [
            Stage00Inputs(),
            Stage10Landmarks(),
            Stage20ExteriorShell(),
            Stage30RegionAtoms(),
            Stage40EmptySpace(),
            Stage50Clustering(),
            Stage55GridRefine05(),   # <--- new
            Stage60SurfaceNormals(),
            Stage70MeshValidate(),
        ]
    )

    pipeline.run(ctx)
    return ctx

```

ribctl/lib/npet2/__init__.py
```py

```

ribctl/lib/npet/kdtree_approach.py
```py
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
import os
from typing import Callable, List, Literal, Optional, Tuple, TypeVar
import open3d as o3d
import numpy as np
from scipy.spatial import cKDTree
import sys
from scipy.spatial import cKDTree
import numpy as np
from Bio.PDB.MMCIFParser import FastMMCIFParser
from typing import List, Tuple

from ribctl.lib.schema.types_ribosome import ResidueSummary

data_dir = os.getenv("DATA_DIR")
sys.dont_write_bytecode = True
from Bio.PDB.Entity import Entity
from Bio.PDB.MMCIFParser import FastMMCIFParser
import numpy as np
from Bio.PDB.Chain import Chain
import pyvista as pv
from pathlib import Path
from pprint import pprint
import subprocess
from matplotlib import pyplot as plt
import open3d as o3d
import pyvista as pv
import numpy as np
import plyfile
import warnings
from ribctl import POISSON_RECON_BIN
import numpy as np
from sklearn.cluster import DBSCAN as skDBSCAN
import requests

warnings.filterwarnings("ignore")
import os


def landmark_constriction_site(rcsb_id: str) -> np.ndarray:
    """
    Fetches the constriction site location for a given RCSB PDB ID.

    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")

    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the constriction site
    """
    url = f"http://localhost:8000/loci/constriction_site?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    response.raise_for_status()  # Raises an exception for 4XX/5XX status codes

    data = response.json()
    return np.array(data["location"])


def landmark_ptc(rcsb_id: str) -> np.ndarray:
    """
    Fetches the peptidyl transferase center (PTC) location for a given RCSB PDB ID.

    Args:
        rcsb_id (str): The RCSB PDB identifier (e.g., "4UG0")

    Returns:
        np.ndarray: Array containing the x, y, z coordinates of the PTC site
    """
    url = f"http://localhost:8000/loci/ptc?rcsb_id={rcsb_id}"
    headers = {"Accept": "application/json"}

    response = requests.get(url, headers=headers)
    response.raise_for_status()
    data = response.json()
    return np.array(data["location"])


def DBSCAN_capture(
    ptcloud: np.ndarray,
    eps,
    min_samples,
    metric: str = "euclidean",
):
    u_EPSILON     = eps
    u_MIN_SAMPLES = min_samples
    u_METRIC      = metric

    print(
        "Running DBSCAN on {} points. eps={}, min_samples={}, distance_metric={}".format(
            len(ptcloud), u_EPSILON, u_MIN_SAMPLES, u_METRIC
        )
    )
    db = skDBSCAN(eps=eps, min_samples=min_samples, metric=metric).fit(ptcloud)
    labels = db.labels_

    CLUSTERS_CONTAINER = {}
    for point, label in zip(ptcloud, labels):
        if label not in CLUSTERS_CONTAINER:
            CLUSTERS_CONTAINER[label] = []
        CLUSTERS_CONTAINER[label].append(point)

    CLUSTERS_CONTAINER = dict(sorted(CLUSTERS_CONTAINER.items()))
    return db, CLUSTERS_CONTAINER


def DBSCAN_pick_largest_cluster(
    clusters_container: dict[int, list], pick_manually: bool = False
) -> tuple[np.ndarray, int]:
    DBSCAN_CLUSTER_ID = 0
    # if pick_manually:
    #     print("-------------------------------")
    #     print("Running Manual Cluster Selection")
    #     picked_id =  int(input("Enter Cluster ID to proceed the reconstruction with\n (options:[{}]):".format(list( clusters_container.keys() ))))
    #     print("Choise cluster # {}".format(picked_id))
    #     if picked_id == -2:
    #         # if picked -2 ==> return largest
    #         for k, v in clusters_container.items():
    #             if int(k) == -1:
    #                 continue
    #             elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
    #                 DBSCAN_CLUSTER_ID = int(k)
    #         return np.array(clusters_container[DBSCAN_CLUSTER_ID])

    #     return np.array(clusters_container[picked_id])

    for k, v in clusters_container.items():
        # print(f"Cluster {k} has {len(v)} points.")
        if int(k) == -1:
            continue
        elif len(v) > len(clusters_container[DBSCAN_CLUSTER_ID]):
            DBSCAN_CLUSTER_ID = int(k)
    return np.array(clusters_container[DBSCAN_CLUSTER_ID]), DBSCAN_CLUSTER_ID


def apply_poisson_reconstruction(
    surf_estimated_ptcloud_path: str,
    output_path: Path,
    recon_depth: int = 6,
    recon_pt_weight: int = 3,
):
    # The documentation can be found at https://www.cs.jhu.edu/~misha/Code/PoissonRecon/Version16.04/ in "PoissonRecon" binary
    print(
        "Rolling Poisson Reconstruction: {} -> {}".format(
            surf_estimated_ptcloud_path, output_path
        )
    )
    command = [
        POISSON_RECON_BIN,
        "--in",
        surf_estimated_ptcloud_path,
        "--out",
        output_path,
        "--depth",
        str(recon_depth),
        "--pointWeight",
        str(recon_pt_weight),
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode == 0:
        data = plyfile.PlyData.read(output_path)
        data.text = True
        ascii_duplicate = output_path.as_posix().split(".")[0] + "_ascii.ply"
        data.write(ascii_duplicate)
        print(">>Wrote {} and the _ascii version.".format(output_path))
    else:
        print(">>Error:", process.stderr)


def ptcloud_convex_hull_points(
    pointcloud: np.ndarray, ALPHA: float, TOLERANCE: float, OFFSET: float
) -> np.ndarray:
    assert pointcloud is not None
    cloud = pv.PolyData(pointcloud)
    grid = cloud.delaunay_3d(
        alpha=ALPHA, tol=TOLERANCE, offset=OFFSET, progress_bar=True
    )
    convex_hull = grid.extract_surface().cast_to_pointset()
    return convex_hull.points


def estimate_normals(
    convex_hull_surface_pts: np.ndarray,
    kdtree_radius=None,
    kdtree_max_nn=None,
    correction_tangent_planes_n=None,
):
    pcd = o3d.geometry.PointCloud()
    pcd.points = o3d.utility.Vector3dVector(convex_hull_surface_pts)
    pcd.estimate_normals(
        search_param=o3d.geometry.KDTreeSearchParamHybrid(
            radius=kdtree_radius, max_nn=kdtree_max_nn
        )
    )
    pcd.orient_normals_consistent_tangent_plane(k=correction_tangent_planes_n)
    # o3d.visualization.draw_geometries([pcd], point_show_normal=True)
    return pcd


T = TypeVar("T")


def generate_voxel_centers(radius: float, height: float, voxel_size: float) -> tuple:
    nx = ny = int(2 * radius / voxel_size) + 1
    nz = int(height / voxel_size) + 1
    x = np.linspace(-radius, radius, nx)
    y = np.linspace(-radius, radius, ny)
    z = np.linspace(0, height, nz)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    voxel_centers = np.column_stack((X.ravel(), Y.ravel(), Z.ravel()))
    return voxel_centers, (X.shape, x, y, z)


def create_point_cloud_mask(
    points: np.ndarray,
    radius: float,
    height: float,
    voxel_size: float = 1.0,
    radius_around_point: float = 2.0,
):
    voxel_centers, (grid_shape, x, y, z) = generate_voxel_centers(
        radius, height, voxel_size
    )
    tree = cKDTree(points)
    indices = tree.query_ball_point(voxel_centers, radius_around_point)

    point_cloud_mask = np.zeros(len(voxel_centers), dtype=bool)
    point_cloud_mask[[i for i, idx in enumerate(indices) if idx]] = True

    point_cloud_mask = point_cloud_mask.reshape(grid_shape)

    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
    cylinder_mask = np.sqrt(X**2 + Y**2) <= radius
    hollow_cylinder = ~cylinder_mask

    final_mask = hollow_cylinder | point_cloud_mask
    return final_mask, (x, y, z)


def get_transformation_to_C0(
    base_point: np.ndarray, axis_point: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute transformation matrices to move arbitrary cylinder to C0 configuration.
    Returns translation vector and rotation matrix.
    """
    # Get cylinder axis vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    # Get rotation that aligns axis_unit with [0, 0, 1]
    z_axis = np.array([0, 0, 1])

    # Use Rodrigues rotation formula to find rotation matrix
    # that rotates axis_unit to z_axis
    if np.allclose(axis_unit, z_axis):
        R = np.eye(3)
    elif np.allclose(axis_unit, -z_axis):
        R = np.diag([1, 1, -1])  # 180-degree rotation around x-axis
    else:
        v = np.cross(axis_unit, z_axis)
        s = np.linalg.norm(v)
        c = np.dot(axis_unit, z_axis)
        v_skew = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        R = np.eye(3) + v_skew + (v_skew @ v_skew) * (1 - c) / (s * s)

    return -base_point, R


def transform_points_to_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_translated = points + translation
    points_transformed = points_translated @ rotation.T

    return points_transformed


def transform_points_from_C0(
    points: np.ndarray, base_point: np.ndarray, axis_point: np.ndarray
) -> np.ndarray:
    translation, rotation = get_transformation_to_C0(base_point, axis_point)
    points_unrotated = points @ rotation
    points_untranslated = points_unrotated - translation

    return points_untranslated


def clip_pointcloud_with_mesh(points: np.ndarray, mesh_path: str) -> np.ndarray:
    """
    Clips a point cloud to keep only points that lie inside a mesh using point-by-point checking.

    Parameters:
        points (np.ndarray): Nx3 array of points to clip
        mesh_path (str): Path to the mesh file (PLY format)

    Returns:
        np.ndarray: Points that lie inside the mesh
    """
    # Load the mesh
    mesh = pv.read(mesh_path)

    # Convert to PolyData if needed
    if not isinstance(mesh, pv.PolyData):
        mesh = mesh.extract_surface()

    # Ensure mesh is triangulated
    if not mesh.is_all_triangles:
        mesh = mesh.triangulate()
    else:
        print("Mesh is triangulated")

    # Initialize mask array
    mask = np.zeros(len(points), dtype=bool)
    print("Got mask", mask.shape)

    # Check each point
    for i, point in enumerate(points):
        # Create a single-point PolyData
        point_data = pv.PolyData(point.reshape(1, 3))
        # Check if point is inside
        selection = mesh.select_enclosed_points(point_data, check_surface=False)
        mask[i] = bool(selection.point_data["SelectedPoints"][0])

        # Progress indicator
        if i % 1000 == 0:
            print(f"Processed {i}/{len(points)} points")

    # Return filtered points
    clipped_points = points[mask]
    print(f"Kept {len(clipped_points)}/{len(points)} points")

    return clipped_points


def verify_mesh_quality(mesh) -> dict:
    """
    Verifies the quality of the input mesh and returns diagnostics.
    """
    stats = {
        "n_points": mesh.n_points,
        "n_faces": mesh.n_faces,
        "is_manifold": mesh.is_manifold,
        "bounds": mesh.bounds,
        "open_edges": mesh.n_open_edges,
    }

    try:
        stats["volume"] = mesh.volume
    except:
        stats["volume"] = None
        print("Warning: Could not compute mesh volume")

    return stats


def clip_pcd_via_ashape(
    pcd: np.ndarray, mesh: pv.PolyData
) -> Tuple[np.ndarray, np.ndarray]:
    points_poly = pv.PolyData(pcd)
    select = points_poly.select_enclosed_points(mesh)
    mask = select["SelectedPoints"]
    ashape_interior = pcd[mask == 1]
    ashape_exterior = pcd[mask == 0]
    return ashape_interior, ashape_exterior


def ribosome_entities(
    rcsb_id: str,
    cifpath: str,
    level=Literal["R"] | Literal["A"],
    skip_nascent_chain: List[str] = [],
) -> list[Entity]:
    structure = FastMMCIFParser(QUIET=True).get_structure(rcsb_id, cifpath)
    residues = []
    for chain in structure.child_list[0]:
        if chain.id in skip_nascent_chain:
            print("Skipping nascent chain ", chain.id)
            continue
        residues.extend(chain)

    for residue in residues:
        if not ResidueSummary.is_canonical(residue.get_resname()):
            residues.remove(residue)

    if level == "R":
        return residues
    elif level == "A":
        atoms = []
        [atoms.extend(a) for a in residues]
        return atoms
    else:
        raise


def is_point_in_cylinder(
    point: np.ndarray,
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
) -> bool:
    point = np.asarray(point)
    base_point = np.asarray(base_point)
    axis_point = np.asarray(axis_point)

    # Calculate cylinder axis direction vector
    axis = axis_point - base_point
    axis_length = np.linalg.norm(axis)
    axis_unit = axis / axis_length

    # Calculate vector from base to point
    point_vector = point - base_point

    # Project vector onto cylinder axis
    projection = np.dot(point_vector, axis_unit)

    # Calculate perpendicular vector from axis to point
    projection_point = base_point + projection * axis_unit
    radial_vector = point - projection_point

    # Calculate radial distance
    radial_distance = np.linalg.norm(radial_vector)

    # Check if point is inside cylinder
    return (radial_distance <= radius) and (-10 <= projection <= height)


def make_cylinder_predicate(
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    position_getter: Callable[[T], np.ndarray],
) -> Callable[[T], bool]:
    def predicate(obj: T) -> bool:
        position = position_getter(obj)
        return is_point_in_cylinder(position, base_point, axis_point, radius, height)

    return predicate


def get_residue_position(residue):
    return residue.center_of_mass()


def _worker_process_chunk(chunk_data):
    """
    Worker function that processes a chunk of residues.
    Takes a tuple of (residue_positions, base_point, axis_point, radius, height, indices)
    Returns indices of residues that are inside the cylinder.
    """
    positions, base_point, axis_point, radius, height, indices = chunk_data

    results = []
    for i, pos in enumerate(positions):
        if is_point_in_cylinder(pos, base_point, axis_point, radius, height):
            results.append(indices[i])

    return results


def filter_residues_parallel(
    residues: List[T],
    base_point: np.ndarray,
    axis_point: np.ndarray,
    radius: float,
    height: float,
    chunk_size: Optional[int] = None,
    max_workers: Optional[int] = None,
) -> List[T]:
    """
    Filter residues in parallel using ProcessPoolExecutor.

    Parameters:
    -----------
    residues : List[T]
        List of residue objects to filter
    base_point : np.ndarray
        Center point of cylinder base
    axis_point : np.ndarray
        Point defining cylinder axis direction
    radius : float
        Radius of cylinder
    height : float
        Height of cylinder
    chunk_size : Optional[int]
        Size of chunks to process in parallel. If None, calculated automatically
    max_workers : Optional[int]
        Maximum number of worker processes. If None, uses CPU count

    Returns:
    --------
    List[T]
        Filtered list of residues whose positions lie within the cylinder
    """
    # Set defaults for parallel processing parameters
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()

    if chunk_size is None:
        # Aim for ~4 chunks per worker for better load balancing
        chunk_size = max(1, len(residues) // (max_workers * 4))

    # Pre-compute all positions and create index mapping
    positions = np.array([get_residue_position(r) for r in residues])

    indices = list(range(len(residues)))
    index_chunks = [
        indices[i : i + chunk_size] for i in range(0, len(indices), chunk_size)
    ]

    # Create data chunks for processing
    chunks_data = [
        (positions[idx], base_point, axis_point, radius, height, idx)
        for idx in index_chunks
    ]

    # Process chunks in parallel
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        results = list(executor.map(_worker_process_chunk, chunks_data))

    # Flatten results and get corresponding residues
    filtered_indices = [idx for chunk_result in results for idx in chunk_result]
    return [residues[i] for i in filtered_indices]


def clip_tunnel_by_chain_proximity(
    tunnel_points: np.ndarray,
    rcsb_id: str,
    cif_path: str,
    chain_id: str = "Y2",
    n_start_residues: int = 7,
    start_proximity_threshold: float = 10.0,  # Tighter threshold for start residues
    rest_proximity_threshold: float = 15.0,  # Wider threshold for rest of chain
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Clips a tunnel point cloud using different proximity thresholds for the start
    and rest of the chain.

    Args:
        tunnel_points: np.ndarray
            Nx3 array of tunnel point cloud coordinates
        rcsb_id: str
            PDB/MMCIF ID of the structure
        cif_path: str
            Path to the MMCIF file
        chain_id: str
            Chain ID of the nascent peptide/reference chain
        n_start_residues: int
            Number of N-terminal residues to apply tighter threshold to
        start_proximity_threshold: float
            Maximum distance (Å) for points near start residues
        rest_proximity_threshold: float
            Maximum distance (Å) for points near rest of chain

    Returns:
        Tuple[np.ndarray, np.ndarray]:
            - Filtered point cloud (points within thresholds)
            - Removed points (points outside thresholds)
    """
    # Parse structure and extract specified chain
    parser = FastMMCIFParser(QUIET=True)
    structure = parser.get_structure(rcsb_id, cif_path)

    # Get the first model and find the specified chain
    model = structure[0]
    try:
        chain = model[chain_id]
    except KeyError:
        raise ValueError(f"Chain {chain_id} not found in structure {rcsb_id}")

    # Calculate centers of mass for each residue in the chain
    residue_centers = []
    residue_indices = []  # Track residue numbers
    for residue in chain:
        coords = np.array([atom.get_coord() for atom in residue])
        center = coords.mean(axis=0)
        residue_centers.append(center)
        residue_indices.append(residue.id[1])  # Get residue number

    residue_centers = np.array(residue_centers)
    residue_indices = np.array(residue_indices)

    # Split centers into start and rest based on residue numbers
    start_centers = residue_centers[residue_indices <= n_start_residues]
    rest_centers = residue_centers[residue_indices > n_start_residues]

    # Build KD-trees for both parts
    start_tree = cKDTree(start_centers)
    rest_tree = cKDTree(rest_centers) if len(rest_centers) > 0 else None

    # Find distances from each tunnel point to nearest residue center
    start_distances, _ = start_tree.query(tunnel_points)
    if rest_tree is not None:
        rest_distances, _ = rest_tree.query(tunnel_points)
    else:
        rest_distances = np.full_like(start_distances, np.inf)

    # Point is kept if it's within threshold of either part
    within_start = start_distances <= start_proximity_threshold
    within_rest = rest_distances <= rest_proximity_threshold
    keep_mask = within_start | within_rest

    # Split points into kept and removed
    kept_points = tunnel_points[keep_mask]
    removed_points = tunnel_points[~keep_mask]

    print(f"Kept {len(kept_points)} points:")
    print(
        f"- {np.sum(within_start)} points within {start_proximity_threshold}Å of first {n_start_residues} residues"
    )
    print(
        f"- {np.sum(within_rest)} points within {rest_proximity_threshold}Å of remaining residues"
    )
    print(f"Removed {len(removed_points)} points")

    return kept_points, removed_points

```


Feel free to ask for any other files...