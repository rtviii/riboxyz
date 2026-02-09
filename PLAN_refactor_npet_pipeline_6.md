Alrighty, i have an NPET extraction pipeline in my riboxyz codebase that i want to clean up and make more robust.

Basically the main issue right now is that we have two main stages of meshing (both results of the dbscan over a voxel cloud) and the first one generally succeeds but the second one generally fails. I'm quite content with the general outline of the pipeline but want to make these two phases robust.

The crux of the matter is that i suppose Poisson Reconstruction algorithm that i used previosuly, when there has only one level of voxel grid refinment (1A) worked somewhat (it still failed sometimes and the steps leading up to it sometimes produced very bad normals and meshes), but now that i have incorporated a furhter voxel grid partitioning (to squeeze our more detail and in fact make the DBSCAN clustering more stable) it fails quite frequently and the task of estimating normals became even hairier shall we say. 

So i'd like to ask for your guidance here as to what we can do to make this bulletproof. Perhaps we can switch entirerly to marching cubes -- the problem here and why i haven't yet done it is that i need a smooth mesh for my colleagues to be able to start a simulation with it and MC produces something quite jagged. But if it is more robust in our context im happy to settle for a "smoothing" post-processing step.

In general this codebase grew pretty adhoc and if you see redundant methods or bad data access (like the `region_atom_xyz_all` below for example) -- let's clean it up, remove bloat etc.

Alrighty let me show you the codebase. Currently it "fails" with the following stuff and produces some pretty shitty output (oversmoothed level0 mehs and empty level1 mesh).

```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  p3 test_npet2.py                                                                                                  [npet_refactor]
[npet2] >>> 00_inputs start
[npet2] <<< 00_inputs done in 1.25s
[npet2] >>> 10_landmarks start
[npet2] <<< 10_landmarks done in 6.91s
[npet2] >>> 20_exterior_shell start
[npet2] <<< 20_exterior_shell done in 0.00s
[npet2] >>> 30_region_atoms start
[30_region_atoms] seed_atoms=15,963 occ_atoms=15,963 occ_chains=53
[npet2] <<< 30_region_atoms done in 12.50s
[npet2] >>> 40_empty_space start
[npet2] <<< 40_empty_space done in 2.64s
[npet2] >>> 50_clustering start
[50_clustering] empty_points n=137,869
Running DBSCAN on 137869 points. eps=5.5, min_samples=600, distance_metric=euclidean
[50_clustering] coarse DBSCAN: 3.47s, 4 clusters
Running DBSCAN on 81006 points. eps=3.5, min_samples=175, distance_metric=euclidean
[50_clustering] refine DBSCAN: 1.13s, 1 clusters
[50_clustering] winner: coarse=81,006 → refine=67,745
[50_clustering] generating mesh for level_0...
Computing 3D Triangulation: 100%|███████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████[00:30<00:00]
[50_clustering]   surface extraction: 30.42s, 11,697 points
[50_clustering]   normal estimation: 0.32s
Rolling Poisson Reconstruction: /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/50_clustering/normals_level_0.ply -> /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/50_clustering/mesh_level_0.ply
>>Wrote /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/50_clustering/mesh_level_0.ply and the _ascii version.
[50_clustering] mesh saved: /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/50_clustering/mesh_level_0.ply
[npet2] <<< 50_clustering done in 36.45s
[npet2] >>> 55_grid_refine start
[55_grid_refine] voxel=0.5Å, ROI pad=10.0Å, atom_r=2.0Å
[55_grid_refine] selected 15,941 atoms near ROI (ALL atoms, prevents interference)
[55_grid_refine] boundary points: 125,822
[55_grid_refine] coarse DBSCAN: 0.52s, 1 clusters
[55_grid_refine] refine DBSCAN: 0.38s, 1 clusters
[55_grid_refine] creating selected void mask for voxel fallback...
[55_grid_refine]   selected mask: 726,627 voxels
[55_grid_refine] DBSCAN refined-grid: coarse=0 → refine=0, final_pts=125,571
[55_grid_refine] generating mesh for level_1...
[55_grid_refine]   normal estimation: 4.22s, 125,571 points
[55_grid_refine]   Poisson reconstruction (o3d, depth=5)...
[55_grid_refine]   Poisson done: 0.15s, 3,032 verts, 6,068 faces
[55_grid_refine] mesh cleanup failed for level_1: 'Missing required input: region_atom_xyz_all'
[npet2] <<< 55_grid_refine done in 14.87s
[npet2] >>> 60_surface_normals start
[60_surface_normals] refined_cluster n=125,571 surface_flag=True
[60_surface_normals] using refined points directly as surface_pts (skip Delaunay)
[60_surface_normals] estimate_normals took 4.24s
[npet2] <<< 60_surface_normals done in 4.26s
[npet2] >>> 70_mesh_validate start
Rolling Poisson Reconstruction: /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/60_surface_normals/surface_normals.ply -> /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/70_mesh_validate/npet2_tunnel_mesh.ply
>>Wrote /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161/stage/70_mesh_validate/npet2_tunnel_mesh.ply and the _ascii version.
[70_mesh_validate] poisson mesh stats: {'n_points': 0, 'n_faces': 0, 'open_edges': 0, 'is_manifold': True, 'bounds': [1e+299, -1e+299, 1e+299, -1e+299, 1e+299, -1e+299]}

Mesh properties at stage: unknown
- Is watertight: True
DONE


[70_mesh_validate] copying comparison meshes...
[70_mesh_validate]   copied level_0 mesh (1.0Å grid)
[70_mesh_validate]   level_1 mesh not found (Poisson likely failed)
[npet2] <<< 70_mesh_validate done in 0.20s
run_dir: /Users/rtviii/dev/riboxyz/NPET2/runs/7K00/20260209_171106_af2eaef57f082161
```

```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.npy|*.ply|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale|NPET2' .  e
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
├── npet2_viewer_usage_examples.md
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_0.md
├── PLAN_refactor_npet_pipeline_1.md
├── PLAN_refactor_npet_pipeline_2.md
├── PLAN_refactor_npet_pipeline_3.md
├── PLAN_refactor_npet_pipeline_4.md
├── PLAN_refactor_npet_pipeline_5.md
├── PLAN_refactor_npet_pipeline_6.md
├── q_entity_filtering.md
├── q_logs_and_params.md
├── q_poisson_still_fails.md
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
│   │   │   │   ├── structure_selection.py
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
└── test_npet2.py
e  [error opening dir]

35 directories, 180 files
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
    impl_version: str = "v1.01"  # bump when you change semantics

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
from typing import List, Literal, Optional, Tuple


@dataclass(frozen=True)
class GridLevelConfig:
    name: str
    voxel_size_A: float

    atom_radius_mode: Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0

    occupancy_backend: Literal["legacy_kdtree", "grid_stamp", "edt", "gpu"] = (
        "legacy_kdtree"
    )
    roi_backend: Literal["full_cylinder", "bbox_from_prev", "tube_from_prev"] = (
        "full_cylinder"
    )


@dataclass(frozen=True)
class RunConfig:
    occupancy_chain_mode: Literal["walls_only", "assembly_all"] = "walls_only"
    occupancy_exclude_trna: bool = True
    occupancy_exclude_auth_asym_ids: Tuple[str, ...] = ()
    # === Region definition ===
    cylinder_radius_A: float = 35.0
    cylinder_height_A: float = 120.0

    # === Stage20: Exterior shell (whole ribosome surface) ===
    # Large alpha because the ribosome is roughly convex
    alpha_d3d_alpha       : float = 200
    alpha_d3d_tol         : float = 10
    alpha_d3d_offset      : float = 3
    alpha_kdtree_radius   : float = 40
    alpha_max_nn          : int   = 60
    alpha_tangent_planes_k: int   = 20
    alpha_poisson_depth   : int   = 6
    alpha_poisson_ptweight: int   = 4
    alpha_fill_holes      : float = 2000

    # === Stage40: Grid levels ===
    grid_levels: List[GridLevelConfig] = field(
        default_factory=lambda: [
            GridLevelConfig(
                name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"
            ),
        ]
    )

    # === Tunnel surface extraction (Delaunay on DBSCAN cluster) ===
    # Small alpha because the tunnel is deeply concave
    tunnel_surface_alpha    : float = 2.0
    tunnel_surface_tolerance: float = 1.0
    tunnel_surface_offset   : float = 2.0

    # === Surface normals (shared by Stage50 mesh, Stage60, Stage55) ===
    normals_radius: float = 15
    normals_max_nn: int = 15
    normals_tangent_k: int = 15

    # === Stage50: DBSCAN clustering on level_0 (1.0A grid) ===
    dbscan_level0_coarse_eps_A: float = 5.5
    dbscan_level0_coarse_min_samples: int = 600
    dbscan_level0_refine_eps_A: float = 3.5
    dbscan_level0_refine_min_samples: int = 175

    # Stage50 mesh generation
    mesh_level0_enable: bool = True
    mesh_level0_poisson_depth: int = 6
    mesh_level0_poisson_ptweight: float = 3.0

    # === Stage55: Grid refinement (0.5A ROI pass) ===
    refine_voxel_size_A: float = 0.5
    refine_roi_pad_A: float = 10.0
    refine_atom_radius_A: float = 2

    refine_keep_within_A: float = 6.0
    refine_occ_close_iters: int = 0
    refine_void_open_iters: int = 1
    refine_forbid_roi_boundary: bool = True

    # DBSCAN on refined grid
    dbscan_level1_coarse_eps_A: float = 3.0
    dbscan_level1_coarse_min_samples: int = 30
    dbscan_level1_refine_eps_A: float = 3.0
    dbscan_level1_refine_min_samples: int = 20

    refine_dbscan_max_points: int = 0
    refine_dbscan_seed: int = 0
    refine_dbscan_max_cluster_stats: int = 25

    # Stage55 mesh generation
    mesh_level1_enable: bool = True
    mesh_level1_poisson_depth: int = 5
    mesh_level1_poisson_ptweight: float = 2

    # === Stage70: Final mesh validation ===
    mesh_poisson_depth: int = 5
    mesh_poisson_ptweight: float = 2.0

    voxel_mesh_fill_holes_A: float = 50.0
    voxel_mesh_smooth_iters: int = 10

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
from datetime import datetime
from typing import Any, Dict


def stable_hash_dict(d: Dict[str, Any]) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def compute_run_id(
    *, 
    rcsb_id: str, 
    pipeline_version: str, 
    inputs_fp: Dict[str, str], 
    config_resolved: Dict[str, Any]
) -> str:
    """
    run_id = timestamp_hash
    
    Format: YYYYMMDD_HHMMSS_<hash16>
    This allows chronological sorting while keeping collision resistance.
    """
    blob = {
        "rcsb_id": rcsb_id.upper(),
        "pipeline_version": pipeline_version,
        "inputs": dict(sorted(inputs_fp.items())),
        "config": config_resolved,
    }
    hash_str = stable_hash_dict(blob)[:16]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{timestamp}_{hash_str}"

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

ribctl/lib/npet2/core/structure_selection.py
```py
# ribctl/lib/npet2/core/structure_selection.py
from __future__ import annotations

from typing import Iterable, Set

_TRNA_HINTS = ("trna", "transfer rna")


def _poly_description(poly) -> str:
    return (getattr(poly, "rcsb_pdbx_description", None) or "").strip()


def _poly_nomenclature_names(poly) -> list[str]:
    noms = getattr(poly, "nomenclature", None) or []
    out: list[str] = []
    for x in noms:
        # PolymerClass enum objects often have .name
        out.append(getattr(x, "name", str(x)))
    return out


def looks_like_trna(poly) -> bool:
    """
    Heuristic: profile-provided description/nomenclature.
    (You said you have a better classifier — you can replace this logic later.)
    """
    desc = _poly_description(poly).lower()
    if any(h in desc for h in _TRNA_HINTS):
        return True

    noms = " ".join(n.lower() for n in _poly_nomenclature_names(poly))
    if "trna" in noms:
        return True

    return False


def ribosome_wall_auth_asym_ids(
    profile,
    *,
    exclude_trna: bool = True,
    extra_exclude: Iterable[str] = (),
) -> Set[str]:
    """
    Define the tunnel 'wall' chains.

    Core rule:
      wall = proteins + rRNAs from the ribosome profile

    Benefits:
      - excludes HOH/ions/ligands/nonpolymers automatically (different chains/entities)
      - includes modified residues inside ribosomal polymers automatically
    """
    proteins = getattr(profile, "proteins", None) or []
    rnas = getattr(profile, "rnas", None) or []

    wall = {p.auth_asym_id for p in proteins} | {r.auth_asym_id for r in rnas}

    if exclude_trna:
        # Defensive: in some datasets a tRNA might be misfiled under rnas/other_polymers.
        others = getattr(profile, "other_polymers", None) or []
        all_polys = list(proteins) + list(rnas) + list(others)
        trna_ids = {p.auth_asym_id for p in all_polys if looks_like_trna(p)}
        wall -= trna_ids

    wall -= set(extra_exclude)
    return wall


def intersect_with_first_assembly(ro, chain_ids: Set[str]) -> Set[str]:
    """
    Optional safety: only use chains present in the first assembly.
    """
    try:
        asm = set(ro.first_assembly_auth_asym_ids())
        return chain_ids & asm
    except Exception:
        return chain_ids

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
import json
from pathlib import Path

from typing import Any, Dict, Tuple

import numpy as np
import pyvista as pv
from scipy import ndimage
import time
import open3d as o3d

from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    transform_points_from_C0,
    estimate_normals,
)

from ribctl.lib.npet2.backends.grid_occupancy import (
    GridSpec,
    occupancy_via_edt,
)


def _make_bbox_grid(lo: np.ndarray, hi: np.ndarray, voxel: float) -> GridSpec:
    lo = np.asarray(lo, dtype=np.float32)
    hi = np.asarray(hi, dtype=np.float32)
    voxel = float(voxel)

    span = hi - lo
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

def _orient_normals_from_selected_mask(
    self,
    ctx: StageContext,
    pcd: o3d.geometry.PointCloud,
    stage_dir: Path,
) -> None:
    """
    Replace pcd.normals with normals derived from the signed distance field (SDF)
    of the selected void component mask. This is robust for branched/loopy tunnels.
    """
    import json
    from scipy import ndimage
    import numpy as np

    # Prefer the selected component mask (the thing you actually want a surface for)
    sel_path = stage_dir / "selected_void_component_mask_level_1.npy"
    if sel_path.exists():
        mask = np.load(sel_path).astype(bool)
    else:
        # fallback
        mask = np.load(stage_dir / "void_mask_level_1.npy").astype(bool)

    # Load grid spec
    spec = json.loads((stage_dir / "grid_spec_level_1.json").read_text())
    grid = GridSpec(
        origin=np.asarray(spec["origin"], dtype=np.float32),
        voxel_size=float(spec["voxel_size_A"]),
        shape=tuple(int(x) for x in spec["shape"]),
    )

    # SDF: positive inside mask, negative outside
    # (distance_transform_edt returns distances in voxel units)
    din  = ndimage.distance_transform_edt(mask)
    dout = ndimage.distance_transform_edt(~mask)
    sdf = din - dout

    # Gradient in ijk index coordinates
    gx, gy, gz = np.gradient(sdf)

    # Points are in world; map to C0 then to ijk
    ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
    constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

    P_w = np.asarray(pcd.points, dtype=np.float32)
    P_c0 = transform_points_to_C0(P_w, ptc, constr).astype(np.float32)

    ijk = _points_to_ijk(grid, P_c0)
    valid = _valid_ijk(grid, ijk)

    # Sample gradient (nearest-voxel)
    N_c0 = np.zeros_like(P_c0, dtype=np.float32)
    ii = ijk[valid, 0]
    jj = ijk[valid, 1]
    kk = ijk[valid, 2]
    N_c0[valid, 0] = gx[ii, jj, kk]
    N_c0[valid, 1] = gy[ii, jj, kk]
    N_c0[valid, 2] = gz[ii, jj, kk]

    # We want normals pointing OUT of the void (toward solid), i.e. toward decreasing sdf
    N_c0[valid] *= -1.0

    # Normalize; if any are zero (rare), keep the existing normal for that point
    N_old = np.asarray(pcd.normals, dtype=np.float32)
    nrm = np.linalg.norm(N_c0, axis=1)
    good = valid & (nrm > 1e-6)
    N_c0[good] /= nrm[good, None]

    # Rotate C0-vector normals back to world. We need the rotation matrix.
    # Import this from your kdtree_approach (it exists in your pasted code)
    from ribctl.lib.npet.kdtree_approach import get_transformation_to_C0
    _, R = get_transformation_to_C0(ptc, constr)  # R maps world->C0 via @ R.T
    N_w = N_c0 @ R  # inverse for vectors: v_world = v_c0 @ R

    N_final = N_old.copy()
    N_final[good] = N_w[good]
    pcd.normals = o3d.utility.Vector3dVector(N_final)


class Stage55GridRefine(Stage):
    key = "55_grid_refine"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A"         : float(c.refine_voxel_size_A),
            "roi_pad_A"            : float(c.refine_roi_pad_A),
            "atom_radius_A"        : float(c.refine_atom_radius_A),
            "keep_within_A"        : float(c.refine_keep_within_A),
            "occ_close_iters"      : int(c.refine_occ_close_iters),
            "void_open_iters"      : int(c.refine_void_open_iters),
            "forbid_roi_boundary"  : bool(c.refine_forbid_roi_boundary),
            "coarse_eps_A"         : float(c.dbscan_level1_coarse_eps_A),
            "coarse_min_samples"   : int(c.dbscan_level1_coarse_min_samples),
            "refine_eps_A"         : float(c.dbscan_level1_refine_eps_A),
            "refine_min_samples"   : int(c.dbscan_level1_refine_min_samples),
            "dbscan_max_points"    : int(c.refine_dbscan_max_points),
            "dbscan_seed"          : int(c.refine_dbscan_seed),
            "mesh_enable"          : bool(getattr(c, "mesh_level1_enable", True)),
            "mesh_poisson_depth"   : int(getattr(c, "mesh_level1_poisson_depth", 8)),
            "mesh_poisson_ptweight": float(getattr(c, "mesh_level1_poisson_ptweight", 0.5)),
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

        refined_world = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)
        
        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = str(ctx.require("alpha_shell_path"))
        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        voxel = float(c.refine_voxel_size_A)
        pad = float(c.refine_roi_pad_A)
        atom_r = float(c.refine_atom_radius_A)
        keep_within_A = float(c.refine_keep_within_A)
        occ_close_iters = int(c.refine_occ_close_iters)
        void_open_iters = int(c.refine_void_open_iters)
        forbid_roi_boundary = bool(c.refine_forbid_roi_boundary)

        eps_c = float(c.dbscan_level1_coarse_eps_A)
        ms_c = int(c.dbscan_level1_coarse_min_samples)
        eps_r = float(c.dbscan_level1_refine_eps_A)
        ms_r = int(c.dbscan_level1_refine_min_samples)

        dbscan_max_points = int(c.refine_dbscan_max_points)
        dbscan_seed = int(c.refine_dbscan_seed)

        print(f"[{self.key}] voxel={voxel}Å, ROI pad={pad}Å, atom_r={atom_r}Å")

        refined_c0 = transform_points_to_C0(refined_world, ptc, constr).astype(np.float32)
        lo = refined_c0.min(axis=0) - pad
        hi = refined_c0.max(axis=0) + pad

        roi_obj = {
            "roi_id": "bbox_pad_stage50",
            "frame": "C0",
            "pad_A": float(pad),
            "lo": [float(x) for x in lo.tolist()],
            "hi": [float(x) for x in hi.tolist()],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
            "source": {"stage": "50_clustering", "artifact": "refined_cluster"},
        }
        (stage_dir / "roi_bbox_c0.json").write_text(json.dumps(roi_obj, indent=2))

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr).astype(np.float32)
        lo_sel = lo - atom_r
        hi_sel = hi + atom_r
        m_atoms = np.all((region_c0 >= lo_sel[None, :]) & (region_c0 <= hi_sel[None, :]), axis=1)
        atoms_roi_c0 = region_c0[m_atoms]
        
        print(f"[{self.key}] selected {atoms_roi_c0.shape[0]:,} atoms near ROI (ALL atoms, prevents interference)")

        grid = _make_bbox_grid(lo, hi, voxel)
        
        cyl = self._cylinder_mask_bbox_grid(
            grid,
            radius_A=float(c.cylinder_radius_A),
            zmin_A=0.0,
            zmax_A=float(c.cylinder_height_A),
        )

        occupied = occupancy_via_edt(atoms_roi_c0, grid, atom_radius_A=atom_r)

        if occ_close_iters > 0:
            occupied = ndimage.binary_closing(occupied, iterations=occ_close_iters)

        occupied = occupied | (~cyl)
        np.save(stage_dir / "occupied_mask_level_1.npy", occupied.astype(np.uint8))

        empty_mask = (~occupied) & cyl
        if empty_mask.sum() == 0:
            raise ValueError(f"[{self.key}] no empty voxels in ROI")

        empty_idx = np.argwhere(empty_mask)
        empty_pts_c0 = _voxel_centers_from_indices(grid, empty_idx).astype(np.float32)

        shell_world = pv.read(alpha_shell_path).triangulate()
        shell_c0 = shell_world.copy(deep=True)
        shell_c0.points = transform_points_to_C0(
            np.asarray(shell_world.points, dtype=np.float32), ptc, constr
        )

        sel = pv.PolyData(empty_pts_c0).select_enclosed_points(shell_c0, check_surface=watertight)
        inside_flags = (np.asarray(sel["SelectedPoints"], dtype=np.int8) == 1)

        void_mask = np.zeros_like(empty_mask, dtype=np.bool_)
        inside_idx = empty_idx[inside_flags]
        if inside_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] no empty voxels inside shell")
        
        void_mask[inside_idx[:, 0], inside_idx[:, 1], inside_idx[:, 2]] = True

        if forbid_roi_boundary:
            void_mask[0, :, :] = False
            void_mask[-1, :, :] = False
            void_mask[:, 0, :] = False
            void_mask[:, -1, :] = False
            void_mask[:, :, 0] = False
            void_mask[:, :, -1] = False

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

        if void_open_iters > 0:
            st = ndimage.generate_binary_structure(3, 1)
            void_mask = ndimage.binary_opening(void_mask, structure=st, iterations=void_open_iters)

        np.save(stage_dir / "void_mask_level_1.npy", void_mask.astype(np.uint8))

        st_er = ndimage.generate_binary_structure(3, 1)
        er = ndimage.binary_erosion(void_mask, structure=st_er, iterations=1)
        boundary_mask = void_mask & (~er)
        boundary_idx = np.argwhere(boundary_mask)
        
        if boundary_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] boundary extraction produced 0 voxels")

        boundary_pts_c0 = _voxel_centers_from_indices(grid, boundary_idx).astype(np.float32)
        boundary_pts_w = transform_points_from_C0(boundary_pts_c0, ptc, constr).astype(np.float32)

        print(f"[{self.key}] boundary points: {boundary_pts_w.shape[0]:,}")

        boundary_pts_c0_cap, boundary_pts_w_cap, idx_cap = self._maybe_cap(
            boundary_pts_c0, boundary_pts_w, dbscan_max_points, dbscan_seed
        )

        tree = cKDTree(refined_c0)

        t0 = time.perf_counter()
        db_coarse = DBSCAN(eps=eps_c, min_samples=ms_c, metric="euclidean", n_jobs=-1)
        labels_coarse = db_coarse.fit_predict(boundary_pts_c0_cap).astype(np.int32)
        dt0 = time.perf_counter() - t0

        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        d_coarse = stage_dir / "coarse"
        self._save_dbscan_pass(
            d_coarse,
            boundary_pts_w_cap,
            labels_coarse,
            eps=eps_c,
            min_samples=ms_c,
        )

        stats_coarse = self._cluster_stats(boundary_pts_c0_cap, labels_coarse, tree)
        best_coarse = self._choose_best_label(stats_coarse)
        
        if best_coarse == -1:
            raise ValueError(
                f"[{self.key}] coarse DBSCAN produced no clusters. "
                f"Try eps={eps_c} → {eps_c*1.5:.1f} or min_samples={ms_c} → {ms_c//2}"
            )

        m_best = labels_coarse == best_coarse
        pts_refine_c0 = boundary_pts_c0_cap[m_best]
        pts_refine_w = boundary_pts_w_cap[m_best]

        t1 = time.perf_counter()
        db_refine = DBSCAN(eps=eps_r, min_samples=ms_r, metric="euclidean", n_jobs=-1)
        labels_refine = db_refine.fit_predict(pts_refine_c0).astype(np.int32)
        dt1 = time.perf_counter() - t1

        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        d_refine = stage_dir / "refine"
        self._save_dbscan_pass(
            d_refine,
            pts_refine_w,
            labels_refine,
            eps=eps_r,
            min_samples=ms_r,
        )

        stats_refine = self._cluster_stats(pts_refine_c0, labels_refine, tree)
        best_refine = self._choose_best_label(stats_refine)
        
        if best_refine == -1:
            raise ValueError(
                f"[{self.key}] refine DBSCAN produced no clusters. "
                f"Try eps={eps_r} → {eps_r*1.5:.1f} or min_samples={ms_r} → {ms_r//2}"
            )

        final_mask = labels_refine == best_refine
        final_surface_w = pts_refine_w[final_mask].astype(np.float32)
        
        if final_surface_w.shape[0] == 0:
            raise ValueError(f"[{self.key}] final cluster is empty")

        np.save(stage_dir / "refined_surface_points_level_1.npy", final_surface_w)

        print(f"[{self.key}] creating selected void mask for voxel fallback...")
        
        final_surface_c0 = transform_points_to_C0(final_surface_w, ptc, constr).astype(np.float32)
        final_ijk = _points_to_ijk(grid, final_surface_c0)
        m_valid = _valid_ijk(grid, final_ijk)
        final_ijk = final_ijk[m_valid]
        
        boundary_selected = np.zeros(grid.shape, dtype=np.bool_)
        if final_ijk.shape[0] > 0:
            boundary_selected[final_ijk[:, 0], final_ijk[:, 1], final_ijk[:, 2]] = True
        
        from scipy.ndimage import binary_fill_holes, binary_dilation
        
        boundary_dilated = binary_dilation(boundary_selected, iterations=2)
        selected_mask = binary_fill_holes(boundary_dilated)
        selected_mask = selected_mask & void_mask
        
        selected_mask_path = stage_dir / "selected_void_component_mask_level_1.npy"
        np.save(selected_mask_path, selected_mask.astype(np.uint8))
        
        print(f"[{self.key}]   selected mask: {selected_mask.sum():,} voxels")
        
        ctx.inputs["selected_void_component_mask_level_1_path"] = str(selected_mask_path)
        ctx.inputs["grid_spec_level_1_path"] = str(stage_dir / "grid_spec_level_1.json")

        spec_obj = {
            "frame": "C0",
            "origin": [float(x) for x in grid.origin.tolist()],
            "voxel_size_A": float(grid.voxel_size),
            "shape": [int(x) for x in grid.shape],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
        }
        (stage_dir / "grid_spec_level_1.json").write_text(json.dumps(spec_obj, indent=2))

        diag = {
            "voxel_size_A": voxel,
            "roi_pad_A": pad,
            "atom_radius_A": atom_r,
            "keep_within_A": keep_within_A,
            "boundary_points_total": int(boundary_pts_c0.shape[0]),
            "boundary_points_used_for_dbscan": int(boundary_pts_c0_cap.shape[0]),
            "dbscan_coarse": {
                "eps_A": eps_c,
                "min_samples": ms_c,
                "best_label": int(best_coarse),
                "n_clusters": int(len(stats_coarse)),
                "clusters": stats_coarse[:25],
            },
            "dbscan_refine": {
                "eps_A": eps_r,
                "min_samples": ms_r,
                "best_label": int(best_refine),
                "n_clusters": int(len(stats_refine)),
                "clusters": stats_refine[:25],
            },
            "final_surface_points": int(final_surface_w.shape[0]),
        }
        (stage_dir / "dbscan_diagnostics.json").write_text(json.dumps(diag, indent=2))

        print(
            f"[{self.key}] DBSCAN refined-grid: "
            f"coarse={best_coarse} → refine={best_refine}, final_pts={final_surface_w.shape[0]:,}"
        )

        ctx.inputs["refined_cluster_surface"] = True
        ctx.inputs["refined_cluster"] = final_surface_w
        ctx.artifacts["refined_surface_points_level_1"] = str(stage_dir / "refined_surface_points_level_1.npy")

        if c.mesh_level1_enable:
            self._generate_mesh(ctx, final_surface_w, level_name="level_1")

    def _cylinder_mask_bbox_grid(
        self, grid: GridSpec, radius_A: float, zmin_A: float, zmax_A: float
    ) -> np.ndarray:
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

    def _maybe_cap(
        self, 
        points_c0: np.ndarray, 
        points_w: np.ndarray, 
        cap: int, 
        seed: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = points_c0.shape[0]
        if cap and n > cap:
            rng = np.random.default_rng(seed)
            idx = rng.choice(n, size=int(cap), replace=False)
            idx.sort()
            return points_c0[idx], points_w[idx], idx
        idx = np.arange(n, dtype=np.int64)
        return points_c0, points_w, idx

    def _cluster_stats(
        self, points_c0: np.ndarray, labels: np.ndarray, tree: cKDTree
    ) -> list[dict]:
        stats = []
        for lab in np.unique(labels):
            if lab == -1:
                continue
            m = labels == lab
            pts = points_c0[m]
            if pts.shape[0] == 0:
                continue
            d, _ = tree.query(pts, k=1)
            stats.append({
                "label": int(lab),
                "size": int(pts.shape[0]),
                "median_dist_to_stage50_A": float(np.median(d)),
                "p05_dist_to_stage50_A": float(np.percentile(d, 5)),
                "p95_dist_to_stage50_A": float(np.percentile(d, 95)),
            })
        stats.sort(key=lambda x: (x["median_dist_to_stage50_A"], -x["size"]))
        return stats

    def _choose_best_label(self, stats: list[dict]) -> int:
        return int(stats[0]["label"]) if stats else -1

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        eps: float,
        min_samples: int,
    ) -> None:
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        clusters_from_labels_func = __import__('ribctl.lib.npet2.backends.clustering_io', fromlist=['clusters_from_labels']).clusters_from_labels
        clusters = clusters_from_labels_func(pts, labels)
        for cid, cpts in clusters.items():
            if cid == -1:
                continue
            if cpts.shape[0] > 0:
                np.save(pass_dir / f"cluster_id{cid}.npy", cpts.astype(np.float32))

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        """Generate mesh from refined boundary points using Poisson reconstruction."""
        import time
        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        print(f"[{self.key}] generating mesh for {level_name}...")

        surface_pts = points

        t1 = time.perf_counter()
        try:
            pcd = estimate_normals(
                surface_pts,
                kdtree_radius=c.normals_radius,
                kdtree_max_nn=c.normals_max_nn,
                correction_tangent_planes_n=c.normals_tangent_k,
            )
        except Exception as e:
            print(f"[{self.key}] normal estimation failed for {level_name}: {e}")
            return
        dt1 = time.perf_counter() - t1
        print(f"[{self.key}]   normal estimation: {dt1:.2f}s, {len(pcd.points):,} points")

        normals_path = stage_dir / f"normals_{level_name}.ply"
        try:
            o3d.io.write_point_cloud(str(normals_path), pcd)
        except Exception as e:
            print(f"[{self.key}] failed to write normals PCD for {level_name}: {e}")
            return

        mesh_path = stage_dir / f"mesh_{level_name}.ply"

        depth = c.mesh_level1_poisson_depth
        print(f"[{self.key}]   Poisson reconstruction (o3d, depth={depth})...")

        t2 = time.perf_counter()
        try:
            mesh_o3d, densities = o3d.geometry.TriangleMesh.create_from_point_cloud_poisson(
                pcd, depth=depth, linear_fit=True
            )
            dt2 = time.perf_counter() - t2
            print(f"[{self.key}]   Poisson done: {dt2:.2f}s, {len(mesh_o3d.vertices):,} verts, {len(mesh_o3d.triangles):,} faces")
        except Exception as e:
            print(f"[{self.key}] o3d Poisson failed for {level_name}: {e}")
            return

        densities = np.asarray(densities)
        density_threshold = np.quantile(densities, 0.01)
        vertices_to_remove = densities < density_threshold
        mesh_o3d.remove_vertices_by_mask(vertices_to_remove)

        o3d.io.write_triangle_mesh(str(mesh_path), mesh_o3d)

        try:
            mesh = pv.read(str(mesh_path))
            mesh = mesh.fill_holes(2000.0)
            mesh = mesh.triangulate()
            mesh = mesh.connectivity(largest=True)

            # Clip mesh back to void -- don't let it eat into atoms
            region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
            mesh = self._clip_mesh_to_atoms(mesh, region_xyz, min_clearance_A=1.5)

            mesh.save(str(mesh_path))

            mesh_path_ascii = stage_dir / f"mesh_{level_name}_ascii.ply"
            mesh.save(str(mesh_path_ascii), binary=False)

            print(f"[{self.key}]   cleaned mesh: {mesh.n_points:,} pts, {mesh.n_faces:,} faces, "
                f"open_edges={mesh.n_open_edges}, manifold={mesh.is_manifold}")

            ctx.store.register_file(
                name=f"mesh_{level_name}",
                stage=self.key,
                type=ArtifactType.PLY_MESH,
                path=mesh_path,
                meta={"level": level_name, "method": "o3d_poisson", "depth": depth},
            )
            print(f"[{self.key}] mesh saved: {mesh_path}")

        except Exception as e:
            print(f"[{self.key}] mesh cleanup failed for {level_name}: {e}")
            try:
                mesh_path.unlink(missing_ok=True)
            except:
                pass
            return

    def _clip_mesh_to_atoms(
        self,
        mesh: pv.PolyData,
        atom_xyz: np.ndarray,
        min_clearance_A: float = 1.5,
    ) -> pv.PolyData:
        """
        Push any mesh vertices that are closer than min_clearance_A to any atom
        back along the atom->vertex direction until they are at min_clearance_A.
        
        This prevents Poisson/smoothing overshoot from eating into atom-occupied space.
        """
        from scipy.spatial import cKDTree

        tree = cKDTree(atom_xyz)
        pts = np.asarray(mesh.points, dtype=np.float64)

        dist, idx = tree.query(pts, k=1)
        violating = dist < min_clearance_A
        n_violations = int(violating.sum())

        if n_violations == 0:
            print(f"[{self.key}]   atom clearance check: all vertices OK (min_clearance={min_clearance_A}A)")
            return mesh

        nearest_atom = atom_xyz[idx[violating]]
        direction = pts[violating] - nearest_atom
        norms = np.linalg.norm(direction, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-8)
        direction = direction / norms

        pts[violating] = nearest_atom + direction * min_clearance_A

        result = mesh.copy()
        result.points = pts.astype(np.float32)

        print(f"[{self.key}]   atom clearance check: pushed {n_violations:,} vertices "
            f"(of {pts.shape[0]:,}) to min {min_clearance_A}A clearance")

        return result
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
from ribctl.lib.npet2.core.structure_selection import intersect_with_first_assembly, ribosome_wall_auth_asym_ids
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
    _valid_ijk,
    _voxel_centers_from_indices,
)
def _residues_from_chain_ids(structure, chain_ids: set[str]):
    model = structure[0]
    residues = []
    for cid in chain_ids:
        if cid not in model:
            continue
        chain = model[cid]
        for r in chain.get_residues():
            # Keep everything inside the polymer chain (including modified residues)
            # Biopython will include modified nucleotides/AAs here.
            if len(getattr(r, "child_list", [])) == 0:
                continue
            residues.append(r)
    return residues

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

        stage_dir        = ctx.store.stage_dir(self.key)
        ptcloud_path     = stage_dir / "ribosome_ptcloud.npy"
        surface_pts_path = stage_dir / "alpha_surface_points.npy"
        normals_pcd_path = stage_dir / "alpha_normals.ply"
        mesh_path        = stage_dir / "alpha_shell.ply"
        quality_path     = stage_dir / "alpha_shell_quality.json"

        # point cloud from cif (legacy)
        # first_assembly_chains = ro.first_assembly_auth_asym_ids()
        # ptcloud = cif_to_point_cloud(
        #     str(cifpath), first_assembly_chains, do_atoms=True
        # ).astype(np.float32)

        ro      = ctx.require("ro")
        profile = ctx.require("profile")
        cifpath = Path(ctx.require("mmcif_path"))

        wall = ribosome_wall_auth_asym_ids(
            profile,
            exclude_trna=bool(getattr(ctx.config, "occupancy_exclude_trna", True)),
            extra_exclude=_tunnel_debris_chains(ctx.rcsb_id, ro, profile),
        )
        wall = intersect_with_first_assembly(ro, wall)

        ptcloud = cif_to_point_cloud(str(cifpath), sorted(wall), do_atoms=True)




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

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        # your existing hardcoded exclusions
        skip = _tunnel_debris_chains(ctx.rcsb_id, ro, profile)

        # plus config-specified exclusions
        skip = list(dict.fromkeys(skip + list(getattr(c, "occupancy_exclude_auth_asym_ids", ()))))


        # ---- choose occupancy chains (the important part) ----
        if getattr(c, "occupancy_chain_mode", "walls_only") == "assembly_all":
            occ_chain_ids = set(ro.first_assembly_auth_asym_ids())
        else:
            occ_chain_ids = ribosome_wall_auth_asym_ids(
                profile,
                exclude_trna=bool(getattr(c, "occupancy_exclude_trna", True)),
                extra_exclude=skip,
            )
            occ_chain_ids = intersect_with_first_assembly(ro, occ_chain_ids)

        # ---- choose seed chains (usually same as occupancy; keep option to diverge later) ----
        seed_chain_ids = set(occ_chain_ids)

        # extract residues from structure once
        structure = ro.assets.biopython_structure()

        residues_seed = _residues_from_chain_ids(structure, seed_chain_ids)
        residues_occ  = _residues_from_chain_ids(structure, occ_chain_ids)

        # spatial filter (cylinder ROI)
        residues_seed = filter_residues_parallel(
            residues=residues_seed,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
        )
        residues_occ = filter_residues_parallel(
            residues=residues_occ,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
        )

        seed_points = np.asarray(
            [atom.get_coord() for r in residues_seed for atom in r.child_list],
            dtype=np.float32,
        )
        occ_points = np.asarray(
            [atom.get_coord() for r in residues_occ for atom in r.child_list],
            dtype=np.float32,
        )

        stage_dir = ctx.store.stage_dir(self.key)

        # keep your existing artifact name for seed points
        out_seed = stage_dir / "region_atom_xyz.npy"
        np.save(out_seed, seed_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_seed,
            meta={"n": int(seed_points.shape[0]), "note": "seed atoms (walls-only chains)"},
        )
        ctx.inputs["region_atom_xyz"] = seed_points

        # NEW: occupancy atoms
        out_occ = stage_dir / "region_atom_xyz_occ.npy"
        np.save(out_occ, occ_points)
        ctx.store.register_file(
            name="region_atom_xyz_occ",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_occ,
            meta={"n": int(occ_points.shape[0]), "note": "occupancy atoms (walls-only chains)"},
        )
        ctx.inputs["region_atom_xyz_occ"] = occ_points

        # Useful debug provenance
        (stage_dir / "occupancy_chain_ids.json").write_text(
            __import__("json").dumps(
                {
                    "occupancy_chain_mode": getattr(c, "occupancy_chain_mode", "walls_only"),
                    "exclude_trna": bool(getattr(c, "occupancy_exclude_trna", True)),
                    "skip_chains": skip,
                    "occupancy_chain_ids": sorted(list(occ_chain_ids)),
                    "seed_chain_ids": sorted(list(seed_chain_ids)),
                },
                indent=2,
            )
        )

        print(
            f"[{self.key}] seed_atoms={seed_points.shape[0]:,} occ_atoms={occ_points.shape[0]:,} "
            f"occ_chains={len(occ_chain_ids)}"
        )



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

        from ribctl.lib.npet2.backends.grid_occupancy import (
            make_cylinder_grid,
            cylinder_mask,
            occupancy_via_edt,
            empty_points_from_mask,
            save_grid_npy,
            get_occupied_voxel_centers,
        )

        c = ctx.config
        
        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        
        # Use filtered atoms for clustering seed reference
        region_xyz_filtered = np.asarray(ctx.require("region_atom_xyz"), dtype=np.float32)
        
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr)

        shell = pv.read(alpha_shell_path)
        if not isinstance(shell, pv.PolyData):
            shell = shell.extract_surface()
        if not shell.is_all_triangles:
            shell = shell.triangulate()

        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        stage_dir = ctx.store.stage_dir(self.key)

        clip_note = {
            "alpha_shell_path": str(alpha_shell_path),
            "alpha_shell_watertight": watertight,
            "clipping_mode": "select_enclosed_points(check_surface=True)"
            if watertight
            else "select_enclosed_points(check_surface=False) [fallback; shell not watertight]",
            "occupancy_atoms": "region_atom_xyz_all (includes ALL chains)",
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
                )
                cyl = np.broadcast_to(cyl2d, grid.shape)

                occ = occ | (~cyl)

                empty_mask = ~occ
                empty_c0 = empty_points_from_mask(grid, empty_mask & cyl)

                save_grid_npy(
                    grid, occ, stage_dir / f"occupancy_grid_{gl.name}", compress=False
                )
                save_grid_npy(
                    grid,
                    (empty_mask & cyl),
                    stage_dir / f"empty_mask_{gl.name}",
                    compress=False,
                )

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
                    pass

            else:
                raise ValueError(
                    f"Grid level {gl.name}: unsupported backend {backend} "
                    f"(supported: legacy_kdtree, edt)"
                )

            empty_world = transform_points_from_C0(empty_c0, ptc, constr).astype(
                np.float32
            )

            p_pre = stage_dir / f"empty_points_{gl.name}_preclip.npy"
            np.save(p_pre, empty_world)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}_preclip",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=p_pre,
                meta={"voxel_size_A": gl.voxel_size_A, "n": int(empty_world.shape[0])},
            )

            if empty_world.shape[0] == 0:
                inside = empty_world
            else:
                pts_poly = pv.PolyData(empty_world)
                sel = pts_poly.select_enclosed_points(shell, check_surface=watertight)
                inside = empty_world[sel["SelectedPoints"] == 1].astype(np.float32)

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
    """
    DBSCAN clustering on level_0 (coarse grid, typically 1.0Å).
    
    Two-pass strategy:
      1. Coarse DBSCAN: merge regions, bridge gaps
      2. Refine DBSCAN: tighten on largest cluster from pass 1
    
    Optionally generates a mesh from the refined cluster.
    
    Outputs:
      stage/50_clustering/
        coarse/{points.npy, labels.npy, cluster_*.npy, index.json}
        refine/{points.npy, labels.npy, cluster_*.npy, index.json}
        largest_cluster.npy
        refined_cluster.npy
        mesh_level_0.ply (if enabled)
    """
    key = "50_clustering"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "coarse_eps_A": c.dbscan_level0_coarse_eps_A,
            "coarse_min_samples": c.dbscan_level0_coarse_min_samples,
            "refine_eps_A": c.dbscan_level0_refine_eps_A,
            "refine_min_samples": c.dbscan_level0_refine_min_samples,
            "mesh_enable": bool(getattr(c, "mesh_level0_enable", True)),
            "mesh_poisson_depth": int(getattr(c, "mesh_level0_poisson_depth", 6)),
            "mesh_poisson_ptweight": int(getattr(c, "mesh_level0_poisson_ptweight", 3)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import time

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        # Input: empty points from Stage40
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)
        if empty_pts.ndim != 2 or empty_pts.shape[1] != 3 or empty_pts.shape[0] == 0:
            raise ValueError(f"[{self.key}] empty_points must be (N,3) and non-empty, got {empty_pts.shape}")

        print(f"[{self.key}] empty_points n={empty_pts.shape[0]:,}")

        # -----------------------
        # PASS 1: Coarse DBSCAN
        # -----------------------
        t0 = time.perf_counter()
        db_coarse, clusters_coarse = DBSCAN_capture(
            empty_pts, 
            c.dbscan_level0_coarse_eps_A, 
            c.dbscan_level0_coarse_min_samples
        )
        labels_coarse = np.asarray(db_coarse.labels_, dtype=np.int32)
        dt0 = time.perf_counter() - t0
        
        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        self._save_dbscan_pass(
            stage_dir / "coarse",
            empty_pts,
            labels_coarse,
            clusters_coarse,
            eps=c.dbscan_level0_coarse_eps_A,
            min_samples=c.dbscan_level0_coarse_min_samples,
        )

        # Pick largest cluster from coarse
        largest, largest_id = DBSCAN_pick_largest_cluster(clusters_coarse)
        largest = np.asarray(largest, dtype=np.float32)
        if largest.shape[0] == 0:
            raise ValueError(f"[{self.key}] largest cluster is empty")

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
        # PASS 2: Refine DBSCAN
        # -----------------------
        t1 = time.perf_counter()
        db_refine, clusters_refine = DBSCAN_capture(
            largest,
            c.dbscan_level0_refine_eps_A,
            c.dbscan_level0_refine_min_samples,
        )
        labels_refine = np.asarray(db_refine.labels_, dtype=np.int32)
        dt1 = time.perf_counter() - t1
        
        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        self._save_dbscan_pass(
            stage_dir / "refine",
            largest,
            labels_refine,
            clusters_refine,
            eps=c.dbscan_level0_refine_eps_A,
            min_samples=c.dbscan_level0_refine_min_samples,
        )

        # Pick winner from refine
        refined, refined_id = DBSCAN_pick_largest_cluster(clusters_refine)
        refined = np.asarray(refined, dtype=np.float32)
        if refined.shape[0] == 0:
            raise ValueError(f"[{self.key}] refined cluster is empty")

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0])},
        )

        # Output for downstream
        ctx.inputs["refined_cluster"] = refined
        ctx.inputs["largest_cluster"] = largest

        print(f"[{self.key}] winner: coarse={largest.shape[0]:,} → refine={refined.shape[0]:,}")

        # -----------------------
        # Optional: Mesh level_0
        # -----------------------
        if c.mesh_level0_enable:
            self._generate_mesh(ctx, refined, level_name="level_0")

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        clusters_dict: dict[int, list],
        eps: float,
        min_samples: int,
    ) -> None:
        """Save DBSCAN pass artifacts."""
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        # Index for quick inspection
        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        for lab, plist in clusters_dict.items():
            lab = int(lab)
            if lab == -1:  # skip noise
                continue
            arr = np.asarray(plist, dtype=np.float32)
            if arr.size > 0:
                np.save(pass_dir / f"cluster_id{lab}.npy", arr)


    # In ribctl/lib/npet2/stages/legacy_minimal.py, Stage50Clustering

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        """Generate mesh from tunnel point cloud using Poisson reconstruction."""
        import time
        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        print(f"[{self.key}] generating mesh for {level_name}...")

        # Surface extraction via tight alpha shape (alpha=2, NOT 200)
        t0 = time.perf_counter()
        try:
            surface_pts = ptcloud_convex_hull_points(
                points,
                ALPHA=c.tunnel_surface_alpha,
                TOLERANCE=c.tunnel_surface_tolerance,
                OFFSET=c.tunnel_surface_offset,
            ).astype(np.float32)
        except Exception as e:
            print(f"[{self.key}] surface extraction failed for {level_name}: {e}")
            return
        dt0 = time.perf_counter() - t0
        print(f"[{self.key}]   surface extraction: {dt0:.2f}s, {surface_pts.shape[0]:,} points")

        # Normal estimation
        t1 = time.perf_counter()
        try:
            pcd = estimate_normals(
                surface_pts,
                kdtree_radius=c.normals_radius,
                kdtree_max_nn=c.normals_max_nn,
                correction_tangent_planes_n=c.normals_tangent_k,
            )
        except Exception as e:
            print(f"[{self.key}] normal estimation failed for {level_name}: {e}")
            return
        dt1 = time.perf_counter() - t1
        print(f"[{self.key}]   normal estimation: {dt1:.2f}s")

        # Write normals PCD
        normals_path = stage_dir / f"normals_{level_name}.ply"
        o3d.io.write_point_cloud(str(normals_path), pcd)

        # Poisson reconstruction
        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        try:
            apply_poisson_reconstruction(
                str(normals_path),
                mesh_path,
                recon_depth=c.mesh_level0_poisson_depth,
                recon_pt_weight=c.mesh_level0_poisson_ptweight,
            )
        except Exception as e:
            print(f"[{self.key}] poisson reconstruction failed for {level_name}: {e}")
            return

        if not mesh_path.exists():
            print(f"[{self.key}] poisson did not produce mesh file for {level_name}")
            return

        # Cleanup mesh
        try:
            mesh = pv.read(str(mesh_path))
            mesh = mesh.fill_holes(2000.0)
            mesh = mesh.connectivity(largest=True).triangulate()
            mesh.save(str(mesh_path))
        except Exception as e:
            print(f"[{self.key}] mesh cleanup failed for {level_name}: {e}")
            return

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name},
        )

        print(f"[{self.key}] mesh saved: {mesh_path}")


class Stage60SurfaceNormals(Stage):
    key = "60_surface_normals"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "tunnel_surface_alpha": c.tunnel_surface_alpha,
            "tunnel_surface_tolerance": c.tunnel_surface_tolerance,
            "tunnel_surface_offset": c.tunnel_surface_offset,
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

        if surface_flag:
            surface_pts = refined
            print(
                "[60_surface_normals] using refined points directly as surface_pts (skip Delaunay)"
            )
        else:
            t0 = time.perf_counter()
            surface_pts = ptcloud_convex_hull_points(
                refined, 
                c.tunnel_surface_alpha,       # was c.surface_alpha
                c.tunnel_surface_tolerance,   # was c.surface_tolerance
                c.tunnel_surface_offset,      # was c.surface_offset
            ).astype(np.float32)
            dt = time.perf_counter() - t0
            print(
                f"[60_surface_normals] delaunay_3d+extract_surface took {dt:,.2f}s surface_pts n={surface_pts.shape[0]:,}"
            )

        p_surface = stage_dir / "surface_points.npy"
        np.save(p_surface, surface_pts)
        ctx.store.register_file(
            name="surface_points",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_surface,
            meta={"n": int(surface_pts.shape[0])},
        )

        t1 = time.perf_counter()
        pcd = estimate_normals(
            surface_pts,
            kdtree_radius=c.normals_radius,
            kdtree_max_nn=c.normals_max_nn,
            correction_tangent_planes_n=c.normals_tangent_k,
        )
        dt1 = time.perf_counter() - t1
        print(f"[60_surface_normals] estimate_normals took {dt1:,.2f}s")

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
            "voxel_fill_holes_A": float(getattr(c, "voxel_mesh_fill_holes_A", 50.0)),
            "voxel_smooth_iters": int(getattr(c, "voxel_mesh_smooth_iters", 10)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import numpy as np
        import pyvista as pv
        from scipy.spatial import cKDTree

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        mesh_path = stage_dir / "npet2_tunnel_mesh.ply"

        def _mesh_stats(m: pv.PolyData) -> dict:
            return {
                "n_points": int(m.n_points),
                "n_faces": int(m.n_faces),
                "open_edges": int(m.n_open_edges),
                "is_manifold": bool(m.is_manifold),
                "bounds": [float(x) for x in m.bounds],
            }

        method_used = None
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

        if mesh_path.exists():
            try:
                m = pv.read(str(mesh_path))
                st = _mesh_stats(m)
                print(f"[70_mesh_validate] poisson mesh stats: {st}")
                watertight = validate_mesh_pyvista(m)
                if watertight:
                    method_used = "poisson"
                    
                    # Clean up Poisson mesh too
                    m = m.connectivity(largest=True)
                    m.save(str(mesh_path))
                    
                    mesh_path_ascii = stage_dir / "npet2_tunnel_mesh_ascii.ply"
                    m.save(str(mesh_path_ascii), binary=False)
                    
                    mesh_path_ascii = stage_dir / "npet2_tunnel_mesh_ascii.ply"
                    try:
                        m.save(str(mesh_path_ascii), binary=False)
                    except:
                        pass
                    
                    ctx.store.register_file(
                        name="tunnel_mesh",
                        stage=self.key,
                        type=ArtifactType.PLY_MESH,
                        path=mesh_path,
                        meta={"watertight": True, "method": "poisson"},
                    )
                    ctx.inputs["tunnel_mesh_path"] = str(mesh_path)
                else:
                    print("[70_mesh_validate] poisson mesh not watertight; falling back to voxel meshing")
            except Exception as e:
                print(f"[70_mesh_validate] failed reading/validating poisson mesh; falling back: {e}")
        else:
            print("[70_mesh_validate] poisson did not produce a mesh file; falling back to voxel meshing")

        if method_used != "poisson":
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

            vol = np.load(mask_p).astype(np.float32)
            if vol.ndim != 3:
                raise ValueError(f"[70_mesh_validate] voxel volume must be 3D, got {vol.shape}")

            vol_pad = np.pad(vol, 1, constant_values=0)
            origin_pad = origin - voxel

            img = pv.ImageData(
                dimensions=vol_pad.shape,
                spacing=(voxel, voxel, voxel),
                origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
            )
            img.point_data["void"] = vol_pad.ravel(order="F")

            surf_c0 = img.contour(isosurfaces=[0.5], scalars="void").triangulate()
            if surf_c0.n_points == 0 or surf_c0.n_faces == 0:
                raise ValueError("[70_mesh_validate] voxel contour produced empty surface")

            surf_c0 = surf_c0.clean(tolerance=0.0)

            fill_holes_A = float(getattr(c, "voxel_mesh_fill_holes_A", 50.0))
            try:
                surf_c0 = surf_c0.fill_holes(fill_holes_A)
            except Exception:
                pass

            smooth_iters = int(getattr(c, "voxel_mesh_smooth_iters", 10))
            if smooth_iters > 0:
                try:
                    surf_c0 = surf_c0.smooth(n_iter=smooth_iters)
                except Exception:
                    pass

            # Largest component LAST -- after fill_holes and smooth
            surf_c0 = surf_c0.connectivity(largest=True)

            pts_c0 = np.asarray(surf_c0.points, dtype=np.float32)
            pts_w = transform_points_from_C0(pts_c0, ptc, constr).astype(np.float32)
            surf_w = surf_c0.copy(deep=True)
            surf_w.points = pts_w

            try:
                surf_w = surf_w.compute_normals(
                    auto_orient_normals=True, consistent_normals=True
                )
            except Exception:
                pass
            # Clip mesh to atom clearance (both Poisson and voxel paths)
            try:
                region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
                from scipy.spatial import cKDTree
                tree = cKDTree(region_xyz)
                pts = np.asarray(surf_w.points, dtype=np.float64)
                dist, idx = tree.query(pts, k=1)
                violating = dist < 1.5
                if violating.sum() > 0:
                    nearest = region_xyz[idx[violating]]
                    direction = pts[violating] - nearest
                    norms = np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-8)
                    pts[violating] = nearest + (direction / norms) * 1.5
                    surf_w.points = pts.astype(np.float32)
                    print(f"[{self.key}]   pushed {int(violating.sum()):,} vertices to 1.5A atom clearance")
            except Exception as e:
                print(f"[{self.key}]   atom clearance clip failed: {e}")

            surf_w.save(str(mesh_path))

            mesh_path_ascii = stage_dir / "npet2_tunnel_mesh_ascii.ply"
            try:
                surf_w.save(str(mesh_path_ascii), binary=False)
                print(f"[{self.key}] saved ASCII mesh: {mesh_path_ascii}")
            except Exception:
                try:
                    import plyfile
                    data = plyfile.PlyData.read(str(mesh_path))
                    data.text = True
                    data.write(str(mesh_path_ascii))
                    print(f"[{self.key}] saved ASCII mesh (via plyfile): {mesh_path_ascii}")
                except Exception as e:
                    print(f"[{self.key}] failed to save ASCII mesh: {e}")

            st2 = _mesh_stats(surf_w)
            print(f"[70_mesh_validate] voxel mesh stats: {st2}")

            watertight = validate_mesh_pyvista(surf_w)
            if not watertight:
                raise ValueError("Final mesh is not watertight (voxel fallback also failed)")

            method_used = "voxel_contour"

        ctx.store.register_file(
            name="tunnel_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"watertight": True, "method": method_used},
        )
        ctx.inputs["tunnel_mesh_path"] = str(mesh_path)

        print(f"[{self.key}] copying comparison meshes...")
        import shutil
        
        try:
            stage50_dir = ctx.store.run_dir / "stage" / "50_clustering"
            mesh_l0_src = stage50_dir / "mesh_level_0.ply"
            if mesh_l0_src.exists():
                mesh_l0_dst = stage_dir / "comparison_mesh_level_0.ply"
                shutil.copy2(mesh_l0_src, mesh_l0_dst)
                
                mesh_l0_src_ascii = stage50_dir / "mesh_level_0_ascii.ply"
                if mesh_l0_src_ascii.exists():
                    mesh_l0_dst_ascii = stage_dir / "comparison_mesh_level_0_ascii.ply"
                    shutil.copy2(mesh_l0_src_ascii, mesh_l0_dst_ascii)
                
                ctx.store.register_file(
                    name="comparison_mesh_level_0",
                    stage=self.key,
                    type=ArtifactType.PLY_MESH,
                    path=mesh_l0_dst,
                    meta={"source": "50_clustering", "voxel_size_A": 1.0},
                )
                print(f"[{self.key}]   copied level_0 mesh (1.0Å grid)")
        except Exception as e:
            print(f"[{self.key}]   failed to copy level_0 mesh: {e}")
        
        try:
            stage55_dir = ctx.store.run_dir / "stage" / "55_grid_refine"
            mesh_l1_src = stage55_dir / "mesh_level_1.ply"
            if mesh_l1_src.exists():
                mesh_l1_dst = stage_dir / "comparison_mesh_level_1.ply"
                shutil.copy2(mesh_l1_src, mesh_l1_dst)
                
                mesh_l1_src_ascii = stage55_dir / "mesh_level_1_ascii.ply"
                if mesh_l1_src_ascii.exists():
                    mesh_l1_dst_ascii = stage_dir / "comparison_mesh_level_1_ascii.ply"
                    shutil.copy2(mesh_l1_src_ascii, mesh_l1_dst_ascii)
                
                ctx.store.register_file(
                    name="comparison_mesh_level_1",
                    stage=self.key,
                    type=ArtifactType.PLY_MESH,
                    path=mesh_l1_dst,
                    meta={"source": "55_grid_refine", "voxel_size_A": 0.5},
                )
                print(f"[{self.key}]   copied level_1 mesh (0.5Å grid)")
            else:
                print(f"[{self.key}]   level_1 mesh not found (Poisson likely failed)")
        except Exception as e:
            print(f"[{self.key}]   failed to copy level_1 mesh: {e}")
```

ribctl/lib/npet2/__init__.py
```py

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
from ribctl.lib.npet2.stages.grid_refine import Stage55GridRefine
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
            Stage55GridRefine(),   # <--- new
            Stage60SurfaceNormals(),
            Stage70MeshValidate(),
        ]
    )

    pipeline.run(ctx)
    return ctx

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
    recon_pt_weight: float = 3.0,  # was int
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

