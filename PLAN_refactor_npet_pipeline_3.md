Alrighty, i have an NPET extraction pipeline in my riboxyz codebase that i want to improve/organize slightly (limit its integration with the actual riboxyz application to only a few key points/coordinates) and then heavily refactor and optimize (bubble up the main computational parameters, implement variable grid resolution and more flexbile parametrization of operations where appropriate, gpu operations where possible).

 I already started by migrating the main parts of the codebase to "npet2". It grew somewhat adhoc.. well, very adhoc when i was building it and currently relies on a mix of a bunch of different disjointed packages and a flow of data that's not robust to a few things that happen pretty often (ex. if the exterior mesh of the ribosome is not watertight then we can't use it to clip the interior mesh, which is pretty important). 

There is also a shred of formal reasoning about parameters we pick for various computational tools during this work (alpha shape, poisson reconstruction, dbscan etc.) so to succeed even half the time everything has to be hand-picked "just so" and my current paramset is a result of heavy trial-and-error.

First salvos for me here seem to be to
- optimize the shell/exterior surface production (so i don't have to rerun it every time if it produces a watertight exterior for some parameters to begin with -- seems like we are well set up to do this with the manifest etc.)
- see if we can indeed improve the visualization facilities for coarse/fine/roi grids as well as the clusters dbscan produces (i'd like to save whatever points dbscan finds, all of them, and be able to visualize them in different colors). We should have a mechanism to visualize differnet sets of clusters (since we might run dbscan mutliple times -- both for refinement inside a single grid resolution and also between different grid resolutions).

After we do this we will be well set up to tackle different/adaptive grid resolutions and possible optimizations/speedups (though im happy to hear what your thoughts already as to what those might be...)


Ok let me show you the code now. I should say that `lib/npet` is the old implemnetation (im happy to share more files from there, just ask me ) and `lib/npet2` is the new one that i'm refactoring out. We've made some motion to switch from ckdtree and dbscan to grid occupancy methods to try to speed things up a bit, but the results were more than dissapointing and im not exactly sure why, but my conclusion was that it's too complex to change multiple of these fundamental methods at once. So i'd say let's be careful in this even though we want to improve multiple of these moving pieces...


```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.npy|*.ply|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale' .  e
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
│       ├── 3J7Z
│       │   └── cf6383ef2118e55a
│       │       ├── manifest.json
│       │       └── stage
│       │           ├── 00_inputs
│       │           ├── 10_landmarks
│       │           ├── 20_exterior_shell
│       │           ├── 30_region_atoms
│       │           ├── 40_empty_space
│       │           ├── 50_clustering
│       │           ├── 60_surface_normals
│       │           └── 70_mesh_validate
│       ├── 4UG0
│       │   └── fa715875c60d9bb0
│       │       ├── manifest.json
│       │       └── stage
│       │           ├── 00_inputs
│       │           ├── 10_landmarks
│       │           ├── 20_exterior_shell
│       │           ├── 30_region_atoms
│       │           ├── 40_empty_space
│       │           ├── 50_clustering
│       │           ├── 60_surface_normals
│       │           └── 70_mesh_validate
│       └── 7K00
│           ├── 54fbca8a18715432
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
│           ├── 658c205113590eca
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
│           └── 874f6b953ff20cc2
│               ├── manifest.json
│               └── stage
│                   ├── 00_inputs
│                   ├── 10_landmarks
│                   ├── 20_exterior_shell
│                   ├── 30_region_atoms
│                   ├── 40_empty_space
│                   ├── 50_clustering
│                   ├── 60_surface_normals
│                   └── 70_mesh_validate
├── npet2_viewer_usage_examples.md
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_1.md
├── PLAN_refactor_npet_pipeline_2.md
├── PLAN_refactor_npet_pipeline_3.md
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
│   │   │   │   ├── grid_occupancy.py
│   │   │   │   └── legacy
│   │   │   ├── core
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
└── v_WKwp5u
e  [error opening dir]

100 directories, 178 files
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]                                                                                  [npet_refactor]

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

ribctl/lib/npet2/backends/__init__.py
```py

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

    # refinement plan (kept for later; default only one level)
    grid_levels: List[GridLevelConfig] = field(default_factory=lambda: [
        GridLevelConfig(name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"),
    ])


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

ribctl/lib/npet2/core/pipeline.py
```py
# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import asdict
from typing import Any, Dict, List, Optional

from .types import StageContext


class Stage(ABC):
    """
    A stage consumes/produces via ctx.inputs and ctx.artifacts.
    Stage key is a stable string like "10_landmarks".
    """

    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]:
        ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None:
        ...


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        for stage in self.stages:
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)
            try:
                stage.run(ctx)
                ctx.store.end_stage(stage.key, success=True)
            except Exception as e:
                ctx.store.end_stage(stage.key, success=False, note=str(e))
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

ribctl/lib/npet2/stages/legacy_minimal.py
```py
from __future__ import annotations
from pathlib import Path
from typing import Any, Dict, List, Tuple
import numpy as np
import pyvista as pv
import open3d as o3d

from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

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
        ptcloud = cif_to_point_cloud(str(cifpath), first_assembly_chains, do_atoms=True).astype(np.float32)
        np.save(ptcloud_path, ptcloud)
        ctx.store.register_file(name="ribosome_ptcloud", stage=self.key, type=ArtifactType.NUMPY, path=ptcloud_path)

        # surface points
        surface_pts = quick_surface_points(ptcloud, c.alpha_d3d_alpha, c.alpha_d3d_tol, c.alpha_d3d_offset).astype(np.float32)
        np.save(surface_pts_path, surface_pts)
        ctx.store.register_file(name="alpha_surface_points", stage=self.key, type=ArtifactType.NUMPY, path=surface_pts_path)

        # normal estimation (legacy)
        normal_estimated_pcd = fast_normal_estimation(surface_pts, c.alpha_kdtree_radius, c.alpha_max_nn, c.alpha_tangent_planes_k)

        # robust-ish normal orientation: outward
        center = normal_estimated_pcd.get_center()
        normal_estimated_pcd.orient_normals_towards_camera_location(camera_location=center)
        normal_estimated_pcd.normals = o3d.utility.Vector3dVector(-np.asarray(normal_estimated_pcd.normals))

        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        ctx.store.register_file(name="alpha_normals_pcd", stage=self.key, type=ArtifactType.PLY_PCD, path=normals_pcd_path)

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
        ctx.store.register_file(name="alpha_shell_quality", stage=self.key, type=ArtifactType.JSON, path=quality_path)

        ctx.store.register_file(name="alpha_shell_mesh", stage=self.key, type=ArtifactType.PLY_MESH, path=mesh_path)
        ctx.inputs["alpha_shell_path"] = str(mesh_path)
        ctx.inputs["alpha_shell_watertight"] = bool(watertight)


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
            max_workers=1,          # <- stops process fan-out
            chunk_size=5000,        # <- optional; irrelevant if max_workers=1
        )


        filtered_points = np.asarray(
            [atom.get_coord() for r in filtered_residues for atom in r.child_list],
            dtype=np.float32,
        )

        out = ctx.store.stage_dir(self.key) / "region_atom_xyz.npy"
        np.save(out, filtered_points)
        ctx.store.register_file(name="region_atom_xyz", stage=self.key, type=ArtifactType.NUMPY, path=out, meta={"n": int(filtered_points.shape[0])})

        ctx.inputs["region_atom_xyz"] = filtered_points


class Stage40EmptySpace(Stage):
    key = "40_empty_space"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "grid_levels": [{"name": gl.name, "voxel_size_A": gl.voxel_size_A, "backend": gl.occupancy_backend} for gl in c.grid_levels],
            "radius_A": c.cylinder_radius_A,
            "height_A": c.cylinder_height_A,
        }

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        region_xyz = np.asarray(ctx.require("region_atom_xyz"), dtype=np.float32)
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        # For now, we only support legacy_kdtree backend (but loop is in place)
        last_empty = None

        for gl in c.grid_levels:
            if gl.occupancy_backend != "legacy_kdtree":
                raise ValueError(f"Grid level {gl.name}: unsupported backend {gl.occupancy_backend} (only legacy_kdtree implemented)")

            # C0 transform
            region_c0 = transform_points_to_C0(region_xyz, ptc, constr)

            # occupancy mask
            mask, (x, y, z) = create_point_cloud_mask(
                region_c0,
                radius=c.cylinder_radius_A,
                height=c.cylinder_height_A,
                voxel_size=gl.voxel_size_A,
                radius_around_point=gl.uniform_atom_radius_A,
            )

            # empty voxel centers (C0)
            idx = np.where(~mask)
            empty_c0 = np.column_stack((x[idx[0]], y[idx[1]], z[idx[2]])).astype(np.float32)

            # back to world
            empty_world = transform_points_from_C0(empty_c0, ptc, constr).astype(np.float32)

            # clip to alpha shell interior
            shell = pv.read(alpha_shell_path)
            sel = pv.PolyData(empty_world).select_enclosed_points(shell)
            inside = empty_world[sel["SelectedPoints"] == 1].astype(np.float32)

            # save artifacts per level
            stage_dir = ctx.store.stage_dir(self.key)
            out = stage_dir / f"empty_points_{gl.name}.npy"
            np.save(out, inside)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={"voxel_size_A": gl.voxel_size_A, "n": int(inside.shape[0])},
            )

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
        }

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)

        # First DBSCAN
        _, clusters = DBSCAN_capture(empty_pts, c.dbscan_eps_A, c.dbscan_min_samples)
        largest, largest_id = DBSCAN_pick_largest_cluster(clusters)

        # Refinement DBSCAN
        _, refined_clusters = DBSCAN_capture(largest, c.refine_eps_A, c.refine_min_samples)
        refined, refined_id = DBSCAN_pick_largest_cluster(refined_clusters)

        stage_dir = ctx.store.stage_dir(self.key)

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest)
        ctx.store.register_file(name="largest_cluster", stage=self.key, type=ArtifactType.NUMPY, path=p_largest, meta={"cluster_id": int(largest_id), "n": int(largest.shape[0])})

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(name="refined_cluster", stage=self.key, type=ArtifactType.NUMPY, path=p_refined, meta={"cluster_id": int(refined_id), "n": int(refined.shape[0])})

        ctx.inputs["refined_cluster"] = refined


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
        c = ctx.config
        refined = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)

        # surface points
        surface_pts = ptcloud_convex_hull_points(refined, c.surface_alpha, c.surface_tolerance, c.surface_offset).astype(np.float32)

        stage_dir = ctx.store.stage_dir(self.key)
        p_surface = stage_dir / "surface_points.npy"
        np.save(p_surface, surface_pts)
        ctx.store.register_file(name="surface_points", stage=self.key, type=ArtifactType.NUMPY, path=p_surface, meta={"n": int(surface_pts.shape[0])})

        # normals
        pcd = estimate_normals(
            surface_pts,
            kdtree_radius=c.normals_radius,
            kdtree_max_nn=c.normals_max_nn,
            correction_tangent_planes_n=c.normals_tangent_k,
        )

        p_normals = stage_dir / "surface_normals.ply"
        o3d.io.write_point_cloud(str(p_normals), pcd)
        ctx.store.register_file(name="surface_normals_pcd", stage=self.key, type=ArtifactType.PLY_PCD, path=p_normals)

        ctx.inputs["normals_pcd_path"] = str(p_normals)


class Stage70MeshValidate(Stage):
    key = "70_mesh_validate"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {"poisson_depth": c.mesh_poisson_depth, "poisson_ptweight": c.mesh_poisson_ptweight}

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        normals_pcd_path = ctx.require("normals_pcd_path")

        stage_dir = ctx.store.stage_dir(self.key)
        mesh_path = stage_dir / "npet2_tunnel_mesh.ply"

        apply_poisson_reconstruction(
            str(normals_pcd_path),
            mesh_path,
            recon_depth=c.mesh_poisson_depth,
            recon_pt_weight=c.mesh_poisson_ptweight,
        )

        watertight = validate_mesh_pyvista(mesh_path)
        if not watertight:
            raise ValueError("Final mesh is not watertight")

        ctx.store.register_file(name="tunnel_mesh", stage=self.key, type=ArtifactType.PLY_MESH, path=mesh_path, meta={"watertight": True})
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
from ribctl.lib.npet2.core.settings import NPET2_RUNS_ROOT
from ribctl.lib.npet2.core.store import LocalRunStore
from ribctl.lib.npet2.core.types import StageContext
from ribctl.lib.npet2.core.pipeline import Pipeline

from ribctl.lib.npet2.stages.bootstrap import Stage00Inputs, Stage10Landmarks
from ribctl.lib.npet2.stages.legacy_minimal import (
    Stage20ExteriorShell,
    Stage30RegionAtoms,
    Stage40EmptySpace,
    Stage50Clustering,
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

    pipeline = Pipeline(
        [
            Stage00Inputs(),
            Stage10Landmarks(),
            Stage20ExteriorShell(),
            Stage30RegionAtoms(),
            Stage40EmptySpace(),
            Stage50Clustering(),
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

ribctl/lib/npet/npet_pipeline.py
```py
from pathlib import Path
from typing import Dict, List, Optional, Any

import numpy as np

from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.kdtree_approach import (
    filter_residues_parallel,
    landmark_constriction_site,
    landmark_ptc,
    ribosome_entities,
)
from ribctl.lib.npet.pipeline.setup_stage import SetupStage
from ribctl.lib.npet.pipeline.ptc_identification_stage import PTCIdentificationStage
from ribctl.lib.npet.pipeline.constriction_identification_stage import (
    ConstrictionIdentificationStage,
)
from ribctl.lib.npet.pipeline.exterior_mesh_stage import AlphaShapeStage
from ribctl.lib.npet.pipeline.landmark_identification_stage import (
    LandmarkIdentificationStage,
)
from ribctl.lib.npet.pipeline.entity_filtering_stage import EntityFilteringStage
from ribctl.lib.npet.pipeline.point_cloud_processing_stage import (
    PointCloudProcessingStage,
)
from ribctl.lib.npet.pipeline.clustering_stage import ClusteringStage
from ribctl.lib.npet.pipeline.refinement_stage import RefinementStage
from ribctl.lib.npet.pipeline.surface_extraction_stage import SurfaceExtractionStage
from ribctl.lib.npet.pipeline.normal_estimation_stage import NormalEstimationStage
from ribctl.lib.npet.pipeline.mesh_reconstruction_stage import MeshReconstructionStage
from ribctl.lib.npet.pipeline.validation_stage import ValidationStage
from ribctl.lib.npet.pipeline_status_tracker import (
    NPETProcessingTracker,
)

def generate_final_cluster_pdb(
    refined_cluster: np.ndarray,
    output_path: str,
    rcsb_id: str
) -> str:
    """
    Generate a PDB file for the final refined cluster used for mesh generation.
    
    Args:
        refined_cluster: The final cluster points array
        output_path: Where to write the PDB file
        rcsb_id: PDB ID for header
        
    Returns:
        str: Path to the generated PDB file
    """
    
    with open(output_path, 'w') as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    FINAL TUNNEL CLUSTER                    {rcsb_id}\n")
        pdb_file.write(f"TITLE     FINAL REFINED CLUSTER FOR {rcsb_id}\n")
        pdb_file.write("REMARK    Final cluster used for mesh generation\n")
        
        atom_num = 1
        connect_records = []
        
        print(f"Writing final cluster with {len(refined_cluster)} points")
        
        # Write each point as a CA atom in chain A, residue 1
        for point_idx, (x, y, z) in enumerate(refined_cluster):
            pdb_file.write(
                f"ATOM  {atom_num:5d}  CA  ALA A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
            )
            
            # Create bonds between consecutive atoms
            if point_idx > 0:
                connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")
            
            atom_num += 1
        
        # Write chain terminator
        pdb_file.write(f"TER   {atom_num:5d}      ALA A   1\n")
        
        # Write all CONECT records
        for connect in connect_records:
            pdb_file.write(connect)
            
        # Write end record
        pdb_file.write("END\n")
    
    print(f"Generated final cluster PDB: {output_path}")
    return output_path

def generate_cylinder_residues_json(
    rcsb_id: str,
    cif_path: str,
    ptc_coords: np.ndarray,
    constriction_coords: np.ndarray,
    output_path: str,
    radius: float = 35.0,
    skip_nascent_chain: List[str] = [],
) -> dict:
    """
    Generate cylinder residues JSON for molstar visualization.

    Args:
        rcsb_id: PDB ID
        cif_path: Path to mmCIF file
        ptc_coords: PTC coordinates from pipeline context
        constriction_coords: Constriction site coordinates from pipeline context
        output_path: Where to write the JSON file
        radius: Cylinder radius in Angstroms
        skip_nascent_chain: Chain IDs to skip

    Returns:
        dict: The cylinder residues data structure
    """
    import json

    # Define cylinder between PTC and constriction site
    base_point = ptc_coords
    axis_point = constriction_coords

    # Calculate height as distance between landmarks plus some buffer
    height = np.linalg.norm(axis_point - base_point) + 10.0

    print(
        f"Cylinder: base={base_point}, axis={axis_point}, radius={radius}, height={height}"
    )

    # Load structure and get residues
    residues = ribosome_entities(
        rcsb_id=rcsb_id,
        cifpath=cif_path,
        level="R",  # Residue level
        skip_nascent_chain=skip_nascent_chain,
    )

    print(f"Loaded {len(residues)} residues from structure")

    # Filter residues that fall within cylinder
    cylinder_residues = filter_residues_parallel(
        residues=residues,
        base_point=base_point,
        axis_point=axis_point,
        radius=radius,
        height=120,
    )

    print(f"Found {len(cylinder_residues)} residues within cylinder")

    # Group by chain and extract residue numbers
    cylinder_data = {}
    for residue in cylinder_residues:
        chain_id = residue.get_parent().id  # auth_asym_id
        residue_num = residue.id[1]  # auth_seq_id

        if chain_id not in cylinder_data:
            cylinder_data[chain_id] = []
        cylinder_data[chain_id].append(residue_num)

    # Sort residue numbers within each chain
    for chain_id in cylinder_data:
        cylinder_data[chain_id].sort()

    print(
        f"Cylinder residues by chain: {dict((k, len(v)) for k, v in cylinder_data.items())}"
    )

    # Write to JSON file
    with open(output_path, "w") as f:
        json.dump(cylinder_data, f, indent=2)

    print(f"Wrote cylinder residues to {output_path}")
    return cylinder_data

def create_npet_mesh(
    rcsb_id: str, log_dir: Optional[Path] = None, force: bool = False
) -> NPETProcessingTracker:
    """
    Creates NPET mesh for a given RCSB ID with comprehensive tracking and logging.

    This function orchestrates the execution of all pipeline stages, handling setup,
    error management, and tracking of the overall process.

    Args:
        rcsb_id: The RCSB PDB identifier
        log_dir: Directory to store logs (default: logs)
        force: Whether to force regeneration of existing assets

    Returns:
        NPETProcessingTracker: Processing tracker with detailed status information
    """
    print(f"Creating NPET mesh for {rcsb_id}")
    meshpath = AssetType.NPET_MESH.get_path(rcsb_id)
    if meshpath.exists():
        print(f"NPET mesh already exists for {rcsb_id}")

    rcsb_id = rcsb_id.upper()
    tracker = NPETProcessingTracker(rcsb_id, log_dir)
    
    try:
        # Initialize pipeline context to share data between stages
        context: Dict[str, Any] = {}
        
        # Define stage parameters
        alpha_shape_params = {
            "d3d_alpha": 200,
            "d3d_tol": 10,
            "d3d_offset": 3,
            "kdtree_radius": 40,
            "max_nn": 60,
            "tangent_planes_k": 20,

            "PR_depth": 6,
            "PR_ptweight": 4,
        }
        
        point_cloud_params = {
            "radius": 35,
            "height": 120,
            "voxel_size": 1,
            "atom_size": 2,
        }
        
        clustering_params = {
            "epsilon": 5.5,
            "min_samples": 600,
        }
        
        refinement_params = {
            "epsilon": 3.5,
            "min_samples": 175,
        }
        
        surface_params = {
            "alpha": 2,
            "tolerance": 1,
            "offset": 2,
        }
        
        normal_params = {
            "kdtree_radius": 10,
            "kdtree_max_nn": 15,
            "correction_tangent_planes_n": 10,
        }
        
        mesh_params = {
            "depth": 6,
            "ptweight": 3,
        }
        
        # Run the Setup stage
        setup_stage = SetupStage(rcsb_id, tracker, None, force)
        setup_result = setup_stage.run(context)
        context.update(setup_result)
        
        # Create artifacts directory
        artifacts_dir = Path(context["ro"].assets.paths.dir) / "artifacts"
        artifacts_dir.mkdir(exist_ok=True, parents=True)
        
        # Run the PTC identification stage
        ptc_stage = PTCIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        ptc_result = ptc_stage.run(context)
        context.update(ptc_result)
        
        # Run the constriction identification stage
        constriction_stage = ConstrictionIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        constriction_result = constriction_stage.run(context)
        context.update(constriction_result)
        
        # Run the alpha shape stage
        alpha_stage = AlphaShapeStage(rcsb_id, tracker, artifacts_dir, True, alpha_shape_params)
        alpha_result = alpha_stage.run(context)
        context.update(alpha_result)
        
        # Run the landmark identification stage (verifies PTC and constriction)
        landmark_stage = LandmarkIdentificationStage(rcsb_id, tracker, artifacts_dir, force)
        landmark_result = landmark_stage.run(context)
        context.update(landmark_result)
        
        # Run the entity filtering stage
        entity_filtering_stage = EntityFilteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"]
        )
        entity_result = entity_filtering_stage.run(context)
        context.update(entity_result)
        
        # Run the point cloud processing stage
        point_cloud_stage = PointCloudProcessingStage(
            rcsb_id, tracker, artifacts_dir, force,
            radius=point_cloud_params["radius"],
            height=point_cloud_params["height"],
            voxel_size=point_cloud_params["voxel_size"],
            atom_size=point_cloud_params["atom_size"]
        )
        point_cloud_result = point_cloud_stage.run(context)
        context.update(point_cloud_result)
        
        # Run the clustering stage
        clustering_stage = ClusteringStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=clustering_params["epsilon"],
            min_samples=clustering_params["min_samples"]
        )
        clustering_result = clustering_stage.run(context)
        context.update(clustering_result)
        
        # Run the refinement stage
        refinement_stage = RefinementStage(
            rcsb_id, tracker, artifacts_dir, force,
            epsilon=refinement_params["epsilon"],
            min_samples=refinement_params["min_samples"]
        )
        refinement_result = refinement_stage.run(context)
        context.update(refinement_result)
        
        # Run the surface extraction stage
        surface_stage = SurfaceExtractionStage(
            rcsb_id, tracker, artifacts_dir, force,
            alpha=surface_params["alpha"],
            tolerance=surface_params["tolerance"],
            offset=surface_params["offset"]
        )
        surface_result = surface_stage.run(context)
        context.update(surface_result)
        
        # Run the normal estimation stage
        normal_stage = NormalEstimationStage(
            rcsb_id, tracker, artifacts_dir, force,
            kdtree_radius=normal_params["kdtree_radius"],
            kdtree_max_nn=normal_params["kdtree_max_nn"],
            correction_tangent_planes_n=normal_params["correction_tangent_planes_n"]
        )
        normal_result = normal_stage.run(context)
        context.update(normal_result)
        
        # Run the mesh reconstruction stage
        mesh_stage = MeshReconstructionStage(
            rcsb_id, tracker, artifacts_dir, force,
            depth=mesh_params["depth"],
            ptweight=mesh_params["ptweight"]
        )
        mesh_result = mesh_stage.run(context)
        context.update(mesh_result)
        
        # Run the validation stage
        validation_stage = ValidationStage(rcsb_id, tracker, artifacts_dir, force)
        validation_result = validation_stage.run(context)
        context.update(validation_result)
        
        # Complete pipeline and set final status
        tracker.complete_processing(True, watertight=context.get("watertight", False))
        
    except Exception as e:
        # Handle any uncaught exceptions
        if tracker.current_stage:
            tracker.end_stage(tracker.current_stage, False, error=e)
        tracker.complete_processing(False)
    
    return tracker
```

__scripts/npet2_viz_enhanced.py
```py
# __scripts/npet2_viz_enhanced.py
#!/usr/bin/env python3
"""
Enhanced NPET2 visualization with grid support.

New features:
- Show occupancy grids (occupied vs empty voxels)
- Show grid boundaries
- Show ROI in world coordinates
- Slice viewer mode
"""

from __future__ import annotations

import argparse
import glob
import json
import re
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pyvista as pv


def find_latest_run_dir(runs_root: Path, rcsb: str) -> Path:
    rdir = runs_root / rcsb.upper()
    if not rdir.exists():
        raise FileNotFoundError(f"No runs for {rcsb} under {runs_root}")
    candidates = [p for p in rdir.iterdir() if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No run subdirs for {rcsb} under {rdir}")
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def load_npy(path: Path) -> np.ndarray:
    a = np.load(path)
    a = np.asarray(a, dtype=np.float32)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"{path} is not (N,3): got {a.shape}")
    return a


def load_grid_npy(base_path: Path) -> tuple[dict, np.ndarray]:
    """Load grid data and spec."""
    data_path = base_path.parent / f"{base_path.stem}_data.npy"
    spec_path = base_path.parent / f"{base_path.stem}_spec.json"
    
    data = np.load(data_path)
    spec = json.loads(spec_path.read_text())
    return spec, data


def grid_to_world_points(spec: dict, data: np.ndarray, ptc: np.ndarray, constr: np.ndarray) -> np.ndarray:
    """
    Convert grid voxel centers from C0 to world coords.
    Only returns centers of True voxels.
    """
    from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
    
    indices = np.argwhere(data)
    origin = np.array(spec["origin"], dtype=np.float32)
    voxel_size = float(spec["voxel_size"])
    
    # Voxel centers in C0
    centers_c0 = origin + indices.astype(np.float32) * voxel_size
    
    # Transform to world
    centers_world = transform_points_from_C0(centers_c0, ptc, constr)
    return centers_world


def parse_cluster_id_from_name(path: Path) -> int:
    m = re.search(r"_id(-?\d+)\.npy$", path.name)
    if not m:
        return 0
    return int(m.group(1))


def add_points(
    plotter: pv.Plotter,
    pts: np.ndarray,
    label: str,
    *,
    point_size: int,
    opacity: float = 1.0,
    color: Optional[str] = None,
):
    cloud = pv.PolyData(pts)
    kwargs = {
        "render_points_as_spheres": True,
        "point_size": point_size,
        "opacity": opacity,
        "label": label,
    }
    if color:
        kwargs["color"] = color
    plotter.add_points(cloud, **kwargs)


def add_roi_bbox(plotter: pv.Plotter, roi_json_path: Path, *, label: str = "ROI bbox", color: str = "black"):
    """
    Load ROI from JSON and display in world coordinates.
    """
    from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
    
    roi = json.loads(roi_json_path.read_text())
    
    # Get transform info
    ptc = np.array(roi["transform"]["ptc"], dtype=np.float32)
    constr = np.array(roi["transform"]["constriction"], dtype=np.float32)
    
    # ROI corners in C0
    lo = np.array(roi["lo"], dtype=np.float32)
    hi = np.array(roi["hi"], dtype=np.float32)
    
    # Generate 8 corners of bbox in C0
    corners_c0 = np.array([
        [lo[0], lo[1], lo[2]],
        [hi[0], lo[1], lo[2]],
        [hi[0], hi[1], lo[2]],
        [lo[0], hi[1], lo[2]],
        [lo[0], lo[1], hi[2]],
        [hi[0], lo[1], hi[2]],
        [hi[0], hi[1], hi[2]],
        [lo[0], hi[1], hi[2]],
    ], dtype=np.float32)
    
    # Transform to world
    corners_world = transform_points_from_C0(corners_c0, ptc, constr)
    
    # Create wireframe box using line segments
    # Define edges of a box (12 edges)
    edges = [
        [0, 1], [1, 2], [2, 3], [3, 0],  # bottom face
        [4, 5], [5, 6], [6, 7], [7, 4],  # top face
        [0, 4], [1, 5], [2, 6], [3, 7],  # vertical edges
    ]
    
    for edge in edges:
        line = pv.Line(corners_world[edge[0]], corners_world[edge[1]])
        plotter.add_mesh(line, color=color, line_width=3, label=label if edge == edges[0] else None)


def add_grid_boundary(plotter: pv.Plotter, spec: dict, ptc: np.ndarray, constr: np.ndarray, 
                     label: str = "Grid boundary", color: str = "gray"):
    """
    Draw wireframe showing the full extent of a voxel grid in world coords.
    """
    from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
    
    origin = np.array(spec["origin"], dtype=np.float32)
    shape = spec["shape"]
    voxel_size = float(spec["voxel_size"])
    
    # Grid extent in C0
    extent = origin + np.array(shape, dtype=np.float32) * voxel_size
    
    # 8 corners
    corners_c0 = np.array([
        [origin[0], origin[1], origin[2]],
        [extent[0], origin[1], origin[2]],
        [extent[0], extent[1], origin[2]],
        [origin[0], extent[1], origin[2]],
        [origin[0], origin[1], extent[2]],
        [extent[0], origin[1], extent[2]],
        [extent[0], extent[1], extent[2]],
        [origin[0], extent[1], extent[2]],
    ], dtype=np.float32)
    
    corners_world = transform_points_from_C0(corners_c0, ptc, constr)
    
    edges = [
        [0, 1], [1, 2], [2, 3], [3, 0],
        [4, 5], [5, 6], [6, 7], [7, 4],
        [0, 4], [1, 5], [2, 6], [3, 7],
    ]
    
    for edge in edges:
        line = pv.Line(corners_world[edge[0]], corners_world[edge[1]])
        plotter.add_mesh(line, color=color, line_width=2, opacity=0.5, 
                        label=label if edge == edges[0] else None)


def add_mesh(plotter: pv.Plotter, mesh_path: Path, label: str, *, opacity: float = 0.15):
    mesh = pv.read(str(mesh_path))
    plotter.add_mesh(mesh, opacity=opacity, label=label)


def add_cluster_group(
    plotter: pv.Plotter,
    paths: List[Path],
    group_label: str,
    *,
    point_size: int,
    downsample: int,
    opacity: float,
):
    if not paths:
        return

    all_pts: List[np.ndarray] = []
    all_lbl: List[np.ndarray] = []

    cids = [parse_cluster_id_from_name(p) for p in paths]
    unique_cids = sorted(set(cids))
    cid_to_idx = {cid: i for i, cid in enumerate(unique_cids)}

    for p in paths:
        cid = parse_cluster_id_from_name(p)
        pts = load_npy(p)
        if downsample > 1:
            pts = pts[::downsample]
        all_pts.append(pts)
        all_lbl.append(np.full((pts.shape[0],), cid_to_idx[cid], dtype=np.int32))

    P = np.vstack(all_pts) if all_pts else np.zeros((0, 3), dtype=np.float32)
    L = np.concatenate(all_lbl) if all_lbl else np.zeros((0,), dtype=np.int32)

    poly = pv.PolyData(P)
    poly["cluster"] = L

    plotter.add_points(
        poly,
        scalars="cluster",
        cmap="tab20",
        render_points_as_spheres=True,
        point_size=point_size,
        opacity=opacity,
        label=group_label,
    )


def main():
    ap = argparse.ArgumentParser(description="Enhanced NPET2 visualization")
    ap.add_argument("--runs_root", default="/Users/rtviii/dev/riboxyz/NPET2/runs")
    ap.add_argument("--rcsb", required=True)
    ap.add_argument("--run_dir", default=None)

    # Existing features
    ap.add_argument("--shell", action="store_true", help="Show alpha shell")
    ap.add_argument("--tunnel_mesh", action="store_true", help="Show tunnel mesh")
    ap.add_argument("--roi", action="store_true", help="Show ROI bounding box")
    ap.add_argument("--points0", action="store_true", help="Show level_0 empty points")
    ap.add_argument("--points1", action="store_true", help="Show level_1 empty points")
    ap.add_argument("--stage40_coarse", action="store_true")
    ap.add_argument("--stage50_coarse", action="store_true")
    ap.add_argument("--stage50_refine", action="store_true")
    ap.add_argument("--stage50_components", action="store_true", 
                   help="Show Stage50 connected components")

    ap.add_argument("--largest", action="store_true")
    ap.add_argument("--refined", action="store_true")

    # **NEW: Grid visualization features**
    ap.add_argument("--occupancy0", action="store_true", help="Show level_0 occupied voxels")
    ap.add_argument("--occupancy1", action="store_true", help="Show level_1 occupied voxels")
    ap.add_argument("--grid_bounds0", action="store_true", help="Show level_0 grid boundaries")
    ap.add_argument("--grid_bounds1", action="store_true", help="Show level_1 grid boundaries")
    ap.add_argument("--empty_mask0", action="store_true", help="Show level_0 empty mask voxels")
    ap.add_argument("--empty_mask1", action="store_true", help="Show level_1 empty mask voxels")


    ap.add_argument("--point_size", type=int, default=3)
    ap.add_argument("--downsample", type=int, default=1)
    ap.add_argument("--opacity", type=float, default=1.0)

    args = ap.parse_args()

    runs_root = Path(args.runs_root)
    if args.run_dir:
        run_dir = Path(args.run_dir)
    else:
        run_dir = find_latest_run_dir(runs_root, args.rcsb)

    print(f"Using run: {run_dir}")

    stage = run_dir / "stage"
    st10 = stage / "10_landmarks"
    st20 = stage / "20_exterior_shell"
    st40 = stage / "40_empty_space"
    st50 = stage / "50_clustering"
    st70 = stage / "70_mesh_validate"

    # Load landmarks for coordinate transform
    ptc_json = st10 / "ptc.json"
    constr_json = st10 / "constriction_site.json"
    
    if not ptc_json.exists() or not constr_json.exists():
        raise FileNotFoundError("Landmarks not found; pipeline may not have run")
    
    ptc = np.array(json.loads(ptc_json.read_text())["location"], dtype=np.float32)
    constr = np.array(json.loads(constr_json.read_text())["location"], dtype=np.float32)

    pl = pv.Plotter()

    # Context meshes
    if args.shell:
        alpha_shell = st20 / "alpha_shell.ply"
        if alpha_shell.exists():
            add_mesh(pl, alpha_shell, "alpha_shell", opacity=0.12)

    if args.tunnel_mesh:
        tunnel = st70 / "npet2_tunnel_mesh.ply"
        if tunnel.exists():
            add_mesh(pl, tunnel, "tunnel_mesh", opacity=0.35)

    # ROI
    if args.roi:
        roi_path = st40 / "roi_bbox_c0.json"
        if roi_path.exists():
            add_roi_bbox(pl, roi_path, label="ROI bbox", color="red")

    # **NEW: Grid boundaries**
    if args.grid_bounds0:
        spec_path = st40 / "occupancy_grid_level_0_spec.json"
        if spec_path.exists():
            spec = json.loads(spec_path.read_text())
            add_grid_boundary(pl, spec, ptc, constr, label="Grid L0 bounds", color="cyan")

    if args.grid_bounds1:
        spec_path = st40 / "occupancy_grid_level_1_spec.json"
        if spec_path.exists():
            spec = json.loads(spec_path.read_text())
            add_grid_boundary(pl, spec, ptc, constr, label="Grid L1 bounds", color="magenta")

    # **NEW: Occupancy visualization**
    if args.occupancy0:
        p = st40 / "occupied_voxels_level_0.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "Occupied L0", point_size=args.point_size, opacity=0.3, color="red")

    if args.occupancy1:
        # Level 1 doesn't have pre-computed occupied voxels; compute on the fly
        spec_path = st40 / "occupancy_grid_level_1_spec.json"
        data_path = st40 / "occupancy_grid_level_1_data.npy"
        if spec_path.exists() and data_path.exists():
            spec, data = load_grid_npy(st40 / "occupancy_grid_level_1")
            pts = grid_to_world_points(spec, data, ptc, constr)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "Occupied L1", point_size=args.point_size, opacity=0.3, color="orange")

    # **NEW: Empty mask visualization**
    if args.empty_mask0:
        spec_path = st40 / "empty_mask_level_0_spec.json"
        data_path = st40 / "empty_mask_level_0_data.npy"
        if spec_path.exists() and data_path.exists():
            spec, data = load_grid_npy(st40 / "empty_mask_level_0")
            pts = grid_to_world_points(spec, data, ptc, constr)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "Empty L0", point_size=args.point_size, opacity=args.opacity, color="cyan")

    if args.empty_mask1:
        spec_path = st40 / "empty_mask_level_1_spec.json"
        data_path = st40 / "empty_mask_level_1_data.npy"
        if spec_path.exists() and data_path.exists():
            spec, data = load_grid_npy(st40 / "empty_mask_level_1")
            pts = grid_to_world_points(spec, data, ptc, constr)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "Empty L1", point_size=args.point_size, opacity=args.opacity, color="blue")

    # Existing empty points
    if args.points0:
        p = st40 / "empty_points_level_0.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "empty_points_L0", point_size=args.point_size, opacity=args.opacity, color="lightblue")

    if args.points1:
        p = st40 / "empty_points_level_1.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "empty_points_L1", point_size=args.point_size, opacity=args.opacity, color="darkblue")

    # Cluster groups
    if args.stage40_coarse:
        paths = sorted(Path(p).resolve() for p in glob.glob(str(st40 / "coarse_cluster_*_id*.npy")))
        add_cluster_group(pl, paths, "Stage40 coarse", point_size=args.point_size, 
                         downsample=args.downsample, opacity=args.opacity)

    if args.stage50_coarse:
        paths = sorted(Path(p).resolve() for p in glob.glob(str(st50 / "coarse_cluster_*_id*.npy")))
        add_cluster_group(pl, paths, "Stage50 coarse", point_size=args.point_size,
                         downsample=args.downsample, opacity=args.opacity)

    if args.stage50_refine:
        paths = sorted(Path(p).resolve() for p in glob.glob(str(st50 / "refine_cluster_*_id*.npy")))
        add_cluster_group(pl, paths, "Stage50 refine", point_size=args.point_size,
                         downsample=args.downsample, opacity=args.opacity)



    if args.stage50_components:
        paths = sorted(Path(p).resolve() for p in glob.glob(str(st50 / "component_*_label*.npy")))
        add_cluster_group(pl, paths, "Stage50 components", point_size=args.point_size,
                         downsample=args.downsample, opacity=args.opacity)
    # Winners
    if args.largest:
        p = st50 / "largest_cluster.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "largest_cluster", point_size=args.point_size + 1, opacity=1.0)

    if args.refined:
        p = st50 / "refined_cluster.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[::args.downsample]
            add_points(pl, pts, "refined_cluster", point_size=args.point_size + 2, opacity=1.0)

    pl.add_axes()
    pl.show_grid()
    pl.add_legend(bcolor="white", size=(0.2, 0.2))
    pl.show()


if __name__ == "__main__":
    main()
```

