Alrighty, i have an NPET extraction pipeline in my riboxyz codebase that i want to improve/organize slightly (limit its integration with the actual riboxyz application to only a few key points/coordinates) and then heavily refactor and optimize (bubble up the main computational parameters, implement variable grid resolution and more flexbile parametrization of operations where appropriate, gpu operations where possible).

For that thought i first want to reorganize this fucking mess. I already started by migrating the main parts of the codebase to "npet2". It grew somewhat adhoc.. well, very adhoc when i was building it and currently relies on a mix of a bunch of different disjointed packages and a flow of data that's not robust to a few things that happen pretty often (ex. if the exterior mesh of the ribosome is not watertight then we can't use it to clip the interior mesh, which is pretty important). 

There is also a shred of formal reasoning about parameters we pick for various computational tools during this work (alpha shape, poisson reconstruction, dbscan etc.) so to succeed even half the time everything has to be hand-picked "just so" and my current paramset is a result of heavy trial-and-error.


I also basically want when we clean this shit up a little what i mean to do this:
- to make things like PTC and Constriction site (coordinates vital to our construction) be ingestable either from the local files or from the riboxyz api
- add a configurable atom radius representation -- right now all atoms are uniformly represented with 2A radius sphere
- add an option to create the grid at other resolution than 1 angstrom wide (i'm not even sure if that's what we use currently...). But basically my hunch is that if the voxel grid is made finer we can actually be more robust to various extraneous clusters and have better separation, have better resolutin in the end. Actually this would be a huge and central improvement but will also be incredibly costly for any finer grid so we should reason very carefully and intentianally about this. This would also of course influence our choice of dbscan parameters and i'd like to have at least some notion of formal/back-of-the-envelope correlation between grid size and the dbscan/poisson reconstruction parameters that i could pass on to the users.
- with all of the above in mind -- i'd like to enable this thing to run on gpu and also parallellize whatever can be parallellized... I'm happy to try to move away from the kdtree to octree or whatever, use other langugaes than python even depending on the tradeoffs..


Here is what we did so far with chatgpt as well(summary and a bit of a gameplan,wbut feel free to suggest other directions..):
```
Below is a “handoff summary” you can paste into the next chat. I’m writing it as a structured recap of **what was reorganized**, **what’s working now**, and **what we agreed to do next** for performance and correctness.

---

## NPET2 refactor summary (legacy → staged pipeline with solid run manifests)

### High-level goal

Create **NPET2** as an operational pipeline with:

* **stable run IDs + per-run artifact directories**
* explicit **stage boundaries + interfaces**
* ability to **swap implementations progressively** (legacy KDTree/Poisson/DBSCAN → grid/EDT/GPU etc.)
* fewer floating files: artifacts live under run dirs, and stages own their outputs.

---

## What we reorganized (structure + interfaces)

### 1) Runs, manifests, and artifact storage

**Implemented: `LocalRunStore` + `manifest.json` per run**

* Run directory layout:

  * `NPET2/runs/<RCSB>/<run_id>/manifest.json`
  * `NPET2/runs/<RCSB>/<run_id>/stage/<stage_key>/...`

* `LocalRunStore` responsibilities:

  * create run dirs + stage dirs
  * write / update `manifest.json`
  * register artifacts consistently (`ArtifactRef`: name/type/path/stage/meta/depends_on)
  * stage lifecycle tracking:

    * `begin_stage(stage, params)`
    * `end_stage(stage, success, note)`
    * `finalize(success, error)`

* Important bug fixed earlier:

  * `_rel()` was missing at one point; later implemented in `store.py` so manifest stores relative paths under `run_dir` whenever possible.

* Decision: **NPET2 root is hardcoded** (for now) rather than mixing with existing ribosome asset manager:

  * `"/Users/rtviii/dev/riboxyz/NPET2"`

This gives you reproducibility, inspectability, and safety: you can re-run and compare runs, and every intermediate artifact is captured.

---

### 2) Pipeline execution harness

**Pipeline object runs ordered Stages, logs params/status.**

* `Stage` interface:

  * `key: str` like `"40_empty_space"`
  * `params(ctx) -> dict` (recorded in manifest)
  * `run(ctx)` (writes artifacts, sets `ctx.inputs` for downstream)

* `Pipeline.run(ctx)`:

  * `store.begin_stage(...)`
  * `stage.run(ctx)`
  * `store.end_stage(...)`
  * on exception: record failure in manifest + re-raise

This makes stage replacements safe: swap stage implementations without rewriting orchestration.

---

## What is still “legacy” vs what was replaced

We kept the pipeline operational by leaving most legacy geometric steps in place, but we replaced the slowest “empty space / cavity sampling” part with the new grid method.

### Stage-by-stage status

#### Stage 00 inputs (operational)

* Stores base inputs like `atom_xyz.npy` etc.

#### Stage 10 landmarks (operational)

* Produces: `ptc.json`, `constriction_site.json`

#### Stage 20 exterior shell (legacy retained)

* Legacy pipeline: point cloud → surface points → normals (Open3D) → Poisson reconstruction → mesh fix (PyVista).
* Artifacts:

  * `ribosome_ptcloud.npy`
  * `alpha_surface_points.npy`
  * `alpha_normals.ply`
  * `alpha_shell.ply` + `_ascii`
  * `alpha_shell_quality.json`
* This is still heavy / dependency-heavy and likely a future replacement candidate.

#### Stage 30 region atoms (legacy retained, but stabilized)

* Extract atoms near the tunnel cylinder region.
* Key fix: stop the “million processes” issue by setting `max_workers=1` in legacy parallel residue filtering.
* Artifact:

  * `region_atom_xyz.npy`

#### Stage 40 empty space (major change: replaced KDTree-ish approach with grid+EDT)

* This was the main “performance + correctness” redesign:

  * Instead of computing “empty points” by brute-ish nearest-neighbor checks / kdtrees repeatedly, we compute an **occupancy grid** and use a **distance transform**.
* Output artifacts:

  * `empty_points_level_0.npy` (coarse voxel centers that are empty)
  * `empty_points_level_1.npy` (finer grid computed inside ROI)
  * `roi_bbox_c0.json` (ROI bounding box in C0 frame)
  * plus saved cluster debugging files (optional):

    * `coarse_cluster_XX_idY.npy` (to visualize DBSCAN output)

#### Stage 50 clustering (legacy retained but reorganized for observability)

* Still uses DBSCAN currently, but organizationally improved:

  * Saves top-K clusters for both passes:

    * coarse pass clusters (on empty points)
    * refined pass clusters (on largest cluster)
  * Still produces the “pipeline outputs”:

    * `largest_cluster.npy`
    * `refined_cluster.npy`
* This fixes the “why do I only see 2 clusters?” issue — previously only winner clusters were saved.

#### Stage 60 surface normals (legacy retained)

* Computes surface points from refined cluster and estimates normals.

#### Stage 70 mesh validate (legacy retained)

* Poisson reconstruction again (Open3D), validate watertight mesh.

---

## The key algorithmic change: occupancy grid + EDT (binary grid + distance transform)

### What we do now (Stage 40)

1. Transform region atoms into **C0 coordinate frame** (PTC/constriction aligned).
2. Build a **cylindrical grid** (voxel centers within bounding box / cylinder).
3. Compute **occupancy**:

   * mark voxels as “occupied” if within `atom_radius_A` of any atom.
4. Take complement inside cylinder:

   * `empty_mask = inside_cylinder & (~occupied)`
5. Convert empty voxels → points (voxel centers) → transform back to world → clip to `alpha_shell` interior (PyVista `select_enclosed_points`).

### Why it’s different from what you had

* Old approach: lots of repeated nearest neighbor queries and geometry operations that can behave like (O(N \log N)) at best but tends to blow up depending on sampling density / repeated passes; also lots of Python overhead.
* New approach:

  * convert problem to a regular grid, then solve “distance-to-atoms” via **Euclidean Distance Transform (EDT)** or equivalent grid-distance logic.
  * EDT-style methods are linear-ish in voxel count (and extremely optimized in SciPy / can be GPU’d).

### Why you didn’t see huge speedup yet

* Your runtime is still dominated by:

  * **Poisson recon + normal estimation + mesh ops** (Stages 20/60/70),
  * **DBSCAN** on large point sets (Stage 50),
  * plus PyVista clipping.
* Stage40 can be faster, but total runtime may still look similar until the other big hitters are replaced or reduced (or you tune voxel sizes).

---

## ROI + progressive refinement (agreed strategy)

### The plan we implemented

* **Level 0**: coarse voxel size over full cylinder region.
* Run a cheap clustering (currently DBSCAN) to estimate where the tunnel void likely is.
* Compute **ROI bounding box** around the largest void component in C0:

  * save as `roi_bbox_c0.json`
  * include margin (e.g. 10Å)
* **Level 1**: finer voxel size only inside that ROI (clamped to cylinder bounds).
* Output becomes `empty_points_level_1.npy`, which Stage50 uses.

### Why this matters

* This avoids paying fine-grid voxel cost over the whole cylinder volume.
* Also keeps memory bounded.
* This is the bridge toward octrees / adaptive refinement later, without needing octree immediately.

---

## Current parameterization and what expanded

### Added / clarified parameters

* Cylinder geometry:

  * `cylinder_radius_A` (you suggested ~35Å)
  * `cylinder_height_A` (~120Å)
* Grid refinement levels:

  * `grid_levels = [{name, voxel_size_A, uniform_atom_radius_A}, ...]`
  * Typically:

    * level0: e.g. 2–4Å voxel
    * level1: e.g. 1–2Å voxel (ROI only)
* ROI parameters:

  * margin (currently hardcoded ~10Å; can be config later)
* Cluster debug saving:

  * Stage40 saves top K coarse clusters
  * Stage50 saves top K clusters for coarse + refine passes

---

## Visualization / benchmarking work we added

### Visualization goal

Be able to inspect:

* empty points coarse vs refined
* ROI bbox correctness
* DBSCAN cluster decomposition (not just largest/refined)
* shell mesh / final mesh overlay

### What was added

* New viewer script approach:

  * load multiple `.npy` clouds
  * render shell mesh semi-transparent
  * draw ROI bbox wireframe
  * for clusters: concatenate all cluster points and attach scalar `cluster_id` so PyVista colormaps them distinctly (solves “everything same color”).

This made it possible to verify Stage40 ROI and cluster behavior visually.

---

## “Menagerie reduction” plan (agreed replacements / performance path)

This is the forward-looking plan we agreed on: keep pipeline interfaces but swap internals for speed + fewer dependencies.

### Biggest runtime hotspots still in legacy

1. **Poisson reconstruction / normal estimation** (Open3D)
2. **DBSCAN** on huge point sets
3. **PyVista mesh operations / clipping**

### Recommended replacement directions

#### A) Keep the grid/EDT approach and deepen it

* Replace “points → DBSCAN → mesh” with “grid → connected components → surface extraction”.
* On a voxel grid you can:

  * select largest connected empty component via **connected-components labeling** (no DBSCAN)
  * use morphological operations, thinning, skeletonization if needed
  * extract surface via **marching cubes** directly (scikit-image / PyMCubes)

This can remove DBSCAN entirely and also reduce reliance on Poisson.

#### B) Connected components instead of DBSCAN (next easy win)

* On voxel grids, clustering becomes trivial:

  * define adjacency (6/18/26-neighborhood)
  * find largest component
* That replaces DBSCAN’s parameter sensitivity and heavy computation.

#### C) Mesh extraction alternatives (reduce Open3D Poisson dependence)

* Instead of Poisson reconstruction, use:

  * **marching cubes** on binary occupancy / SDF
  * potentially smoothing then marching cubes
* Tools:

  * `skimage.measure.marching_cubes`
  * or VTK’s marching cubes directly (if already depending on PyVista/VTK)
* This reduces Open3D dependency and avoids Poisson instability warnings.

#### D) GPU enablement path

* Grid/EDT is a natural GPU target:

  * CuPy / cupyx.scipy.ndimage distance transforms
  * PyTorch tensor ops (for distance field approximation / occupancy)
* But we decided: **don’t prematurely GPU** — first remove Poisson/DBSCAN bottlenecks.

#### E) Octree discussion (future, not first)

* Octrees can do adaptive refinement like ROI without fixed levels.
* But they add complexity; since we already have level0→ROI→level1, octree is optional later.
* If adopted:

  * you’d store a sparse occupancy / SDF octree
  * extract cavity surface similarly (marching cubes on adaptive structure or dual contouring)
* Agreed: do EDT grid first; revisit octree only if memory/resolution become limiting.

---

## Operational status at handoff

* Pipeline runs end-to-end for multiple structures.
* Stage40 now produces meaningful empty-space point clouds (level0/level1).
* ROI bbox is persisted in stage40 and viewable.
* DBSCAN clusters can be persisted and visualized (top-K saved).
* Remaining issues are mostly:

  * parameter tuning (voxel sizes, atom radius, ROI margin)
  * replacing DBSCAN and Poisson for major speedups
  * cleaning dependencies (Open3D/PyVista/Trimesh sprawl)

---

## Next concrete steps we agreed on (order of attack)

1. **Replace Stage50 DBSCAN with grid connected-components** (using Stage40 occupancy / emptiness masks, not just sampled points).
2. **Replace mesh construction (Stages 60–70) with marching cubes** on an SDF / occupancy grid for the refined component.
3. Reduce dependence on Open3D by removing Poisson reconstruction path.
4. Only after (1–3), consider:

   * GPU EDT / SDF
   * octree / adaptive sparse grids
   * performance tuning and batch benchmarking.

---

If you want, I can also write a short “design contract” for each stage (inputs/outputs/expected coordinate frame + artifact names) so the next LLM can safely rewrite internals without breaking your manifest + run structure.

```


Ok here is the repo and _some_ of the code...


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
│   ├── npet2_view.py
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
├── f_1UZWoj
├── f_3i8cI4
├── f_6B5VLR
├── f_6WKwJf
├── f_Fsj2S8
├── f_GhB3T2
├── f_kLsi6b
├── f_LXu1Kc
├── f_mYQ35k
├── f_YcZ5Tw
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
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_1.md
├── PLAN_refactor_npet_pipeline_2.md
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
└── v_WKwp5u
e  [error opening dir]

70 directories, 179 files
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
    alpha_d3d_alpha       : float = 200
    alpha_d3d_tol         : float = 10
    alpha_d3d_offset      : float = 3
    alpha_kdtree_radius   : float = 40
    alpha_max_nn          : int   = 60
    alpha_tangent_planes_k: int   = 20
    alpha_poisson_depth   : int   = 6
    alpha_poisson_ptweight: int   = 4
    alpha_fill_holes      : float = 2000

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
    grid_levels = [
        GridLevelConfig(name="level_0", voxel_size_A=1.5, occupancy_backend="grid_edt", uniform_atom_radius_A=2.0),
        GridLevelConfig(name="level_1", voxel_size_A=1.0, occupancy_backend="grid_edt", uniform_atom_radius_A=2.0),
    ]



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

import open3d as o3d
# keep using your transforms + DBSCAN + pyvista shell clip for now
import pyvista as pv
import numpy as np
from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from ribctl.lib.npet2.backends.grid_occupancy import (
    make_cylinder_grid,
    cylinder_mask,
    occupancy_via_edt,
    empty_points_from_mask,
    GridSpec,
)


from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    transform_points_from_C0,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
)


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

    def params(self, ctx):
        c = ctx.config
        return {
            "radius_A": c.cylinder_radius_A,
            "height_A": c.cylinder_height_A,
            "grid_levels": [{"name": gl.name, "voxel_size_A": gl.voxel_size_A} for gl in c.grid_levels],
            "atom_radius_A": c.grid_levels[0].uniform_atom_radius_A if c.grid_levels else 2.0,
        }

    def run(self, ctx):
        c = ctx.config
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        region_xyz = np.asarray(ctx.require("region_atom_xyz"), dtype=np.float32)
        region_c0  = transform_points_to_C0(region_xyz, ptc, constr)

        # We’ll assume at least one level exists
        levels = c.grid_levels
        if not levels:
            raise ValueError("config.grid_levels is empty")

        # Helper: clip points (world) to alpha shell interior
        shell = pv.read(alpha_shell_path)

        def clip_to_shell(world_pts: np.ndarray) -> np.ndarray:
            sel = pv.PolyData(world_pts).select_enclosed_points(shell)
            return world_pts[sel["SelectedPoints"] == 1].astype(np.float32)

        # -------------------------
        # Level 0 (coarse full cyl)
        # -------------------------
        gl0 = levels[0]
        grid0 = make_cylinder_grid(c.cylinder_radius_A, c.cylinder_height_A, gl0.voxel_size_A)
        inside_cyl0 = cylinder_mask(grid0, c.cylinder_radius_A)

        occ0 = occupancy_via_edt(region_c0, grid0, atom_radius_A=gl0.uniform_atom_radius_A)

        # empty in cylinder AND not occupied
        empty0_mask = inside_cyl0 & (~occ0)
        empty0_c0 = empty_points_from_mask(grid0, empty0_mask)

        # back to world + clip to shell
        empty0_world = transform_points_from_C0(empty0_c0, ptc, constr)
        empty0_world = clip_to_shell(empty0_world)

        # save empty0
        stage_dir = ctx.store.stage_dir(self.key)
        p_empty0 = stage_dir / f"empty_points_{gl0.name}.npy"
        np.save(p_empty0, empty0_world)
        ctx.store.register_file(
            name=f"empty_points_{gl0.name}",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_empty0,
            meta={"voxel_A": gl0.voxel_size_A, "n": int(empty0_world.shape[0])},
        )

        # -----------------------------------------
        # Coarse tunnel ROI estimate (cheap DBSCAN)
        # -----------------------------------------
        # NOTE: This is “rough”; it’s just to define ROI for finer grid.
        # You can tune these later or compute from config heuristics.
        eps0 = max(3.0, 5.5 * gl0.voxel_size_A)  # loose-ish scaling
        min0 = 300  # reduce for coarse; you can parameterize later

        _, clusters0 = DBSCAN_capture(empty0_world, eps0, min0)
        coarse_cluster_world, _ = DBSCAN_pick_largest_cluster(clusters0)

        # Save coarse clusters (optionally top K to avoid writing dozens/hundreds)
        items = [(cid, np.asarray(pts, dtype=np.float32)) for cid, pts in clusters0.items()]
        items.sort(key=lambda x: x[1].shape[0], reverse=True)

        K = 10  # change as needed
        for rank, (cid, pts) in enumerate(items[:K]):
            tag = f"coarse_cluster_{rank:02d}_id{cid}"
            p = stage_dir / f"{tag}.npy"
            np.save(p, pts)
            ctx.store.register_file(
                name=tag,
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=p,
                meta={"cluster_id": int(cid), "rank": int(rank), "n": int(pts.shape[0]), "eps": float(eps0), "min_samples": int(min0)},
            )


        # define ROI in C0
        coarse_cluster_c0 = transform_points_to_C0(coarse_cluster_world, ptc, constr)
        lo = coarse_cluster_c0.min(axis=0)
        hi = coarse_cluster_c0.max(axis=0)

        margin = 10.0
        lo -= margin
        hi += margin
        roi = {
            "space": "C0",
            "lo": lo.astype(float).tolist(),
            "hi": hi.astype(float).tolist(),
            "size": (hi - lo).astype(float).tolist(),
            "margin_A": float(margin),
            "level0": {"name": gl0.name, "voxel_A": float(gl0.voxel_size_A)},
        }

        # If you have a level1 configured, log it too (optional)
        if len(levels) > 1:
            roi["level1"] = {"name": levels[1].name, "voxel_A": float(levels[1].voxel_size_A)}

        ctx.store.put_json(name="roi_bbox_c0", stage=self.key, obj=roi)




        # If no next level, we are done: Stage50 will cluster/refine empty0
        if len(levels) == 1:
            ctx.inputs["empty_points"] = empty0_world
            return

        # --------------------------
        # Level 1 (fine inside ROI)
        # --------------------------
        gl1 = levels[1]
        v1 = float(gl1.voxel_size_A)

        # Create a grid spec for ROI bbox intersected with cylinder bounds
        # Keep it axis-aligned in C0 for simplicity.
        # Clamp ROI to cylinder box bounds
        cyl_lo = np.array([-c.cylinder_radius_A, -c.cylinder_radius_A, 0.0], dtype=np.float32)
        cyl_hi = np.array([ c.cylinder_radius_A,  c.cylinder_radius_A, c.cylinder_height_A], dtype=np.float32)
        lo = np.maximum(lo, cyl_lo)
        hi = np.minimum(hi, cyl_hi)

        # Build ROI grid
        nx = int(np.floor((hi[0] - lo[0]) / v1)) + 1
        ny = int(np.floor((hi[1] - lo[1]) / v1)) + 1
        nz = int(np.floor((hi[2] - lo[2]) / v1)) + 1
        grid1 = GridSpec(origin=lo.astype(np.float32), voxel_size=v1, shape=(nx, ny, nz))

        # cylinder mask in ROI grid coordinates
        inside_cyl1 = cylinder_mask(grid1, c.cylinder_radius_A)

        occ1 = occupancy_via_edt(region_c0, grid1, atom_radius_A=gl1.uniform_atom_radius_A)
        empty1_mask = inside_cyl1 & (~occ1)
        empty1_c0 = empty_points_from_mask(grid1, empty1_mask)

        empty1_world = transform_points_from_C0(empty1_c0, ptc, constr)
        empty1_world = clip_to_shell(empty1_world)

        p_empty1 = stage_dir / f"empty_points_{gl1.name}.npy"
        np.save(p_empty1, empty1_world)
        ctx.store.register_file(
            name=f"empty_points_{gl1.name}",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_empty1,
            meta={"voxel_A": gl1.voxel_size_A, "n": int(empty1_world.shape[0]), "roi_margin_A": margin},
        )


        print(f"[Stage40] ROI C0 lo={lo}, hi={hi}, size={hi-lo}")
        print(f"[Stage40] level0 empty={len(empty0_world)}; level1 voxel={gl1.voxel_size_A}A")
        print(f"[Stage40] ROI grid shape={grid1.shape} voxels (~{grid1.shape[0]*grid1.shape[1]*grid1.shape[2]:,})")

        # Output for next stage
        ctx.inputs["empty_points"] = empty1_world



class Stage50Clustering(Stage):
    key = "50_clustering"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "dbscan_eps_A": c.dbscan_eps_A,
            "dbscan_min_samples": c.dbscan_min_samples,
            "refine_eps_A": c.refine_eps_A,
            "refine_min_samples": c.refine_min_samples,
            # organizational knobs:
            "save_top_k": getattr(c, "save_top_k_clusters", 12),
            "save_noise": getattr(c, "save_noise_cluster", False),
        }

    @staticmethod
    def _save_clusters(
        *,
        ctx: StageContext,
        stage_dir: Path,
        clusters: Dict[int, np.ndarray],
        prefix: str,
        eps: float,
        min_samples: int,
        top_k: int,
        save_noise: bool,
    ) -> List[Tuple[int, Path, int]]:
        """
        Save cluster point clouds as .npy so we can inspect DBSCAN behavior.
        Returns [(cluster_id, path, n_points), ...] saved in size-desc order.
        """
        items: List[Tuple[int, np.ndarray]] = []
        for cid, pts in clusters.items():
            if cid == -1 and not save_noise:
                continue
            if pts is None:
                continue
            pts = np.asarray(pts, dtype=np.float32)
            if pts.size == 0:
                continue
            items.append((int(cid), pts))

        # Sort by size desc
        items.sort(key=lambda x: x[1].shape[0], reverse=True)

        saved: List[Tuple[int, Path, int]] = []
        for rank, (cid, pts) in enumerate(items[:top_k]):
            tag = f"{prefix}_cluster_{rank:02d}_id{cid}"
            out = stage_dir / f"{tag}.npy"
            np.save(out, pts)
            ctx.store.register_file(
                name=tag,
                stage=ctx.store.manifest.stages.get(ctx.store.manifest.stages.keys().__iter__().__next__(), None) and "" or "",  # unused
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "cluster_id": int(cid),
                    "rank": int(rank),
                    "n": int(pts.shape[0]),
                    "eps": float(eps),
                    "min_samples": int(min_samples),
                    "prefix": prefix,
                },
            )
            # IMPORTANT: correct stage name for store
            # (register_file requires stage string; set to self.key in caller)
            saved.append((int(cid), out, int(pts.shape[0])))

        return saved

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)

        stage_dir = ctx.store.stage_dir(self.key)

        top_k = getattr(c, "save_top_k_clusters", 12)
        save_noise = getattr(c, "save_noise_cluster", False)

        # -----------------------
        # Pass 1: DBSCAN on empty
        # -----------------------
        _, clusters = DBSCAN_capture(empty_pts, c.dbscan_eps_A, c.dbscan_min_samples)

        # Save top-K clusters for visualization/debug
        items = [(int(cid), np.asarray(pts, dtype=np.float32)) for cid, pts in clusters.items() if pts is not None]
        if not save_noise:
            items = [(cid, pts) for cid, pts in items if cid != -1]
        items = [(cid, pts) for cid, pts in items if pts.size > 0]
        items.sort(key=lambda x: x[1].shape[0], reverse=True)

        for rank, (cid, pts) in enumerate(items[:top_k]):
            tag = f"coarse_cluster_{rank:02d}_id{cid}"
            out = stage_dir / f"{tag}.npy"
            np.save(out, pts)
            ctx.store.register_file(
                name=tag,
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "cluster_id": int(cid),
                    "rank": int(rank),
                    "n": int(pts.shape[0]),
                    "eps": float(c.dbscan_eps_A),
                    "min_samples": int(c.dbscan_min_samples),
                },
            )

        # Pick winner for pipeline continuity
        largest, largest_id = DBSCAN_pick_largest_cluster(clusters)

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest.astype(np.float32))
        ctx.store.register_file(
            name="largest_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_largest,
            meta={"cluster_id": int(largest_id), "n": int(largest.shape[0])},
        )

        # ----------------------------------
        # Pass 2: refinement DBSCAN on winner
        # ----------------------------------
        _, refined_clusters = DBSCAN_capture(largest, c.refine_eps_A, c.refine_min_samples)

        # Save top-K refined clusters too
        ritems = [(int(cid), np.asarray(pts, dtype=np.float32)) for cid, pts in refined_clusters.items() if pts is not None]
        if not save_noise:
            ritems = [(cid, pts) for cid, pts in ritems if cid != -1]
        ritems = [(cid, pts) for cid, pts in ritems if pts.size > 0]
        ritems.sort(key=lambda x: x[1].shape[0], reverse=True)

        for rank, (cid, pts) in enumerate(ritems[:top_k]):
            tag = f"refine_cluster_{rank:02d}_id{cid}"
            out = stage_dir / f"{tag}.npy"
            np.save(out, pts)
            ctx.store.register_file(
                name=tag,
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "cluster_id": int(cid),
                    "rank": int(rank),
                    "n": int(pts.shape[0]),
                    "eps": float(c.refine_eps_A),
                    "min_samples": int(c.refine_min_samples),
                },
            )

        refined, refined_id = DBSCAN_pick_largest_cluster(refined_clusters)

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined.astype(np.float32))
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0])},
        )

        ctx.inputs["refined_cluster"] = refined.astype(np.float32)



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

ribctl/lib/npet/tunnel_asset_manager.py
```py
import os

from ribctl import ASSETS_PATH, RIBETL_DATA, RIBXZ_TEMP_FILES
from ribctl.ribosome_ops import RibosomeOps


class TunnelMeshAssetsManager:

    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()
        RO = RibosomeOps(rcsb_id)
        self.structpath = RO.assets.paths.cif
        pass

    @property
    def cif_struct(self):
        return self.structpath

    @property
    def tunnel_pcd_normal_estimated(self):
        return os.path.join( RIBXZ_TEMP_FILES, "{}_tunnel_pcd_normal_estimated.ply".format(self.rcsb_id) )

    @property
    def tunnel_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_tunnel_half_poisson_recon.ply".format(self.rcsb_id)
        )

    @property
    def ashape_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_half_ashape_watertight.ply".format(self.rcsb_id)
        )

    # These are just to visualize lamps trajectories/pointclouds as pdb files.

    @property
    def lammps_traj_tunnel(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_tunnel_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.xyz".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.xyz".format(self.rcsb_id)
        )
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

ribctl/lib/npet/pipeline_status_tracker.py
```py
import copy
import json
import os
import time
import traceback
from enum import Enum
from pathlib import Path
from typing import Dict, Any, Optional, List

class ProcessingStage(str, Enum):
    """Enum defining the stages of the NPET mesh pipeline"""
    SETUP = "setup"

    PTC_IDENTIFICATION = "ptc_identification"           # New stage
    CONSTRICTION_IDENTIFICATION = "constriction_identification"  # New stage
    ALPHA_SHAPE = "alpha_shape"
    LANDMARK_IDENTIFICATION = "landmark_identification"
    ENTITY_FILTERING = "entity_filtering"
    POINT_CLOUD_PROCESSING = "point_cloud_processing"
    CLUSTERING = "clustering"
    REFINEMENT = "refinement"
    SURFACE_EXTRACTION = "surface_extraction"
    NORMAL_ESTIMATION = "normal_estimation"
    MESH_RECONSTRUCTION = "mesh_reconstruction"
    VALIDATION = "validation"
    COMPLETE = "complete"

class ProcessingStatus(str, Enum):
    """Status of a processing stage"""
    PENDING     = "pending"
    IN_PROGRESS = "in_progress"
    SUCCESS = "success"
    FAILURE = "failure"
    SKIPPED = "skipped"

class NPETProcessingTracker:
    """Tracks the processing status of NPET mesh generation for a structure"""
    
    def __init__(self, rcsb_id: str, output_dir: Optional[Path] = None):
        self.rcsb_id = rcsb_id.upper()
        self.start_time = time.time()
        self.end_time: Optional[float] = None
        self.current_stage: ProcessingStage = ProcessingStage.SETUP
        self.output_dir = output_dir or Path(f"logs")
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize stage tracking
        self.stages: Dict[ProcessingStage, Dict[str, Any]] = {
            stage: {
                "status": ProcessingStatus.PENDING,
                "start_time": None,
                "end_time": None,
                "duration": None,
                "error": None,
                "artifacts": [],
                "parameters": {}  # New field to store parameters
            } for stage in ProcessingStage
        }
        
        # Initialize overall status
        self.status: Dict[str, Any] = {
            "rcsb_id": self.rcsb_id,
            "overall_status": ProcessingStatus.PENDING,
            "duration": None,
            "stages": self.stages,
            "summary": {
                "success": False,
                "failed_stage": None,
                "error_summary": None,
                "watertight": False,
                "artifacts_generated": []
            }
        }
        
        # Save initial status
        self._save_status()
    
    def begin_stage(self, stage: ProcessingStage, parameters: Optional[Dict[str, Any]] = None) -> bool:
        """
        Mark the beginning of a processing stage.
        
        Args:
            stage: The processing stage
            parameters: Optional dictionary of parameters for this stage
            
        Returns:
            bool: True if the stage should be processed, False if it can be skipped
        """
        self.current_stage = stage
        now = time.time()
        
        # Print stage start
        print(f"\n⏳ Starting stage: {stage}")
        
        # Check if we can skip processing based on parameters
        can_skip = False
        if parameters is not None:
            # Get previous parameters and artifacts
            prev_params = self.stages[stage].get("parameters", {})
            prev_artifacts = self.stages[stage].get("artifacts", [])
            
            # Check if parameters match and artifacts exist
            if prev_params and prev_artifacts and prev_params == parameters:
                # Check that artifacts actually exist on disk
                artifacts_exist = all(os.path.exists(artifact) for artifact in prev_artifacts)
                if artifacts_exist:
                    can_skip = True
                    print(f"⏩ Stage {stage} skipped (parameters unchanged, artifacts exist)")
                    
                    # Update status but keep existing data
                    self.stages[stage].update({
                        "status": ProcessingStatus.SUCCESS,
                        "start_time": now,
                        "end_time": now,
                        "duration": 0.0,  # Zero duration for skipped stages
                    })
                    self._save_status()
                    return False  # Skip processing
        
        # Proceed with processing
        self.stages[stage].update({
            "status": ProcessingStatus.IN_PROGRESS,
            "start_time": now,
            "parameters": parameters or {}  # Store new parameters
        })
        self._save_status()
        return True  # Process the stage
    
    def end_stage(self, stage: ProcessingStage, success: bool, 
                artifacts: Optional[List[Path]] = None,
                error: Optional[Exception] = None) -> None:
        """Mark the end of a processing stage"""
        now = time.time()
        start_time = self.stages[stage]["start_time"] or now
        
        status = ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE
        
        error_info = None
        if error:
            error_info = {
                "type": type(error).__name__,
                "message": str(error),
                "traceback": traceback.format_exc()
            }
        
        # Create a dict with the updates, but don't include artifacts unless provided
        update_dict = {
            "status": status,
            "end_time": now,
            "duration": now - start_time,
            "error": error_info,
        }
        
        # Only update the artifacts list if explicitly provided
        if artifacts is not None:
            update_dict["artifacts"] = [str(a) for a in artifacts]
        
        # Apply the updates
        self.stages[stage].update(update_dict)
        
        if not success:
            # Update summary with failure info
            self.status["summary"]["success"] = False
            self.status["summary"]["failed_stage"] = stage
            self.status["summary"]["error_summary"] = str(error) if error else "Unknown error"
            
            # Mark remaining stages as skipped
            all_stages = list(ProcessingStage)
            current_index = all_stages.index(stage)
            for skip_stage in all_stages[current_index+1:]:
                if skip_stage != ProcessingStage.COMPLETE:
                    self.stages[skip_stage]["status"] = ProcessingStatus.SKIPPED
                    
            # Print failure message with emoji
            print(f"❌ Stage FAILED: {stage} - {str(error) if error else 'Unknown error'}")
        else:
            # Print success message with emoji
            duration = now - start_time
            print(f"✅ Stage completed: {stage} ({duration:.2f}s)")
        
        # Print artifact count
        artifact_count = len(self.stages[stage].get("artifacts", []))
        if artifact_count > 0:
            print(f"   {artifact_count} artifacts generated")
        
        self._save_status()
    
    def complete_processing(self, success: bool, watertight: bool = False) -> None:
        """Mark the overall processing as complete"""
        now = time.time()
        self.end_time = now
        total_duration = now - self.start_time
        
        # Get all generated artifacts with debugging
        artifacts = []
        for stage_name, stage_info in self.stages.items():
            stage_artifacts = stage_info.get("artifacts", [])
            if stage_artifacts:
                print(f"Found {len(stage_artifacts)} artifacts in stage {stage_name}: {stage_artifacts}")
                artifacts.extend(stage_artifacts)
        
        print(f"Total artifacts collected: {len(artifacts)}")
        
        self.status.update({
            "overall_status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "duration": total_duration,
            "summary": {
                "success": success,
                "watertight": watertight,
                "artifacts_generated": artifacts
            }
        })
        
        # Update the complete stage
        self.stages[ProcessingStage.COMPLETE].update({
            "status": ProcessingStatus.SUCCESS if success else ProcessingStatus.FAILURE,
            "duration": 0
        })
        
        # Don't print summary here - let the caller handle it
        # self.print_summary()
        
        self._save_status()
    
    def add_artifact(self, stage: ProcessingStage, artifact_path: Path) -> None:
        """Add an artifact to the specified stage"""
        artifact_str = str(artifact_path)
        
        if artifact_str not in self.stages[stage].get("artifacts", []):
            self.stages[stage].setdefault("artifacts", []).append(artifact_str)
            self._save_status()
    
    def print_summary(self) -> None:
        """Print a summary of all stages with their status"""
        print("\n" + "=" * 60)
        print(f"PROCESSING SUMMARY FOR {self.rcsb_id}")
        print("=" * 60)
        
        # Calculate the widest stage name for nice formatting
        stage_width = max(len(str(stage)) for stage in ProcessingStage) + 2
        
        # Print status for each stage
        for stage in ProcessingStage:
            if stage != ProcessingStage.COMPLETE:
                status = self.stages[stage]["status"]
                
                # Select appropriate emoji and color based on status
                emoji = {
                    ProcessingStatus.SUCCESS: "✅",
                    ProcessingStatus.FAILURE: "❌",
                    ProcessingStatus.SKIPPED: "⏩",
                    ProcessingStatus.IN_PROGRESS: "⏳",
                    ProcessingStatus.PENDING: "⏱️",
                }.get(status, "❓")
                
                # Get duration if available
                duration_str = ""
                if self.stages[stage]["duration"] is not None:
                    duration = self.stages[stage]["duration"]
                    if duration < 0.001:  # Skipped stages
                        duration_str = "(skipped)"
                    else:
                        duration_str = f"({duration:.2f}s)"
                
                # Get artifacts if available
                artifact_count = len(self.stages[stage].get("artifacts", []))
                artifact_str = f"{artifact_count} artifacts" if artifact_count > 0 else ""
                
                # Get error if available
                error_str = ""
                if self.stages[stage]["error"]:
                    error_str = f"ERROR: {self.stages[stage]['error']['message']}"
                
                # Format the line
                stage_name = f"{stage}".ljust(stage_width)
                status_name = f"{status.name}".ljust(10)
                duration_str = duration_str.ljust(15)
                artifact_str = artifact_str.ljust(15)
                
                print(f"{emoji} {stage_name} {status_name} {duration_str} {artifact_str} {error_str}")
        
        print("-" * 60)
        
        # Print overall status
        total_duration = self.status.get("duration", 0)
        duration_str = f"Total time: {total_duration:.2f}s" if total_duration else ""
        
        success = self.status["summary"]["success"]
        if success:
            print(f"✅ PIPELINE COMPLETED SUCCESSFULLY  {duration_str}")
            if self.status["summary"]["watertight"]:
                mesh_files = [a for a in self.status["summary"]["artifacts_generated"] 
                            if a.endswith('.ply') and not a.endswith('_normal_estimated_pcd.ply')]
                if mesh_files:
                    print(f"   Watertight mesh: {os.path.basename(mesh_files[0])}")
        else:
            failed_stage = self.status["summary"]["failed_stage"]
            error = self.status["summary"]["error_summary"]
            print(f"❌ PIPELINE FAILED at stage {failed_stage}  {duration_str}")
            if error:
                print(f"   Error: {error}")
        
        print("=" * 60)
    
    def _save_status(self) -> None:
        """Save the current status to a JSON file"""
        log_file = self.output_dir / f"{self.rcsb_id}_processing_log.json"
        
        # Create a simplified version of the status for JSON output
        output_status = copy.deepcopy(self.status)
        
        # Simplify the time representation - remove start_time and end_time
        for stage_name, stage_info in output_status["stages"].items():
            if "start_time" in stage_info:
                del stage_info["start_time"]
            if "end_time" in stage_info:
                del stage_info["end_time"]
        
        with open(log_file, 'w') as f:
            json.dump(output_status, f, indent=2, default=str)
```


Wanna see anything else?