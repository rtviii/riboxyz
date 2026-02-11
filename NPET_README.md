# NPET2 Pipeline

Computes the ribosome exit tunnel (NPET) as a watertight mesh from an mmCIF structure.

## Pipeline stages

### 00_inputs
Loads the structure from mmCIF via the configured StructureProvider. Extracts all atom coordinates and stores them as a numpy array. Computes bounding box stats for downstream stages.

No config parameters.

### 10_landmarks
Computes two reference points that define the tunnel axis:
- **PTC** (peptidyl transferase center) -- the "bottom" of the tunnel
- **Constriction site** -- the narrowest point, roughly mid-tunnel

From these, derives the cylinder z-extents in C0 (canonical orientation where PTC is at z=0 and the axis is along z).

Config parameters:
- `cylinder_height_A` -- total height of the search cylinder along the PTC-Constriction axis. Default 150.
- `cylinder_ptc_extension_A` -- how far below PTC (in angstroms) the cylinder extends. This captures the PTC-tRNA interface region. Default 15.

### 20_exterior_shell
Builds a watertight mesh of the entire ribosome surface using Delaunay 3D + Poisson reconstruction. This shell is later used to clip void space -- only empty voxels *inside* the ribosome are kept.

Results are cached by structure fingerprint, so repeated runs on the same structure skip this stage.

Config parameters:
- `alpha_d3d_alpha` -- Delaunay 3D alpha value. Larger = coarser surface. Default 200.
- `alpha_d3d_tol` -- Delaunay tolerance. Default 10.
- `alpha_d3d_offset` -- Delaunay offset. Default 3.
- `alpha_kdtree_radius` -- KDTree search radius for normal estimation. Default 40.
- `alpha_max_nn` -- max neighbors for normal estimation. Default 60.
- `alpha_tangent_planes_k` -- tangent plane consistency parameter. Default 20.
- `alpha_poisson_depth` -- Poisson reconstruction depth. Higher = finer detail but slower. Default 6.
- `alpha_poisson_ptweight` -- Poisson point weight. Default 4.
- `alpha_fill_holes` -- hole-filling threshold on the shell mesh. Default 2000.

### 30_region_atoms
Selects which atoms define "occupied space" (tunnel walls) within the search cylinder.

Inclusion policy:
- **Included**: ribosomal proteins + rRNAs, including all modified residues (pseudouridine, methylated bases, etc.) since they are covalently part of the wall.
- **Excluded**: water molecules (HOH), ions (Mg, K, etc.), nonpolymer ligands (antibiotics, spermidine, paromomycin), tRNAs, manually specified chains, known debris chains.

Writes `atom_selection_policy.json` with full provenance of what was included/excluded and why.

Config parameters:
- `occupancy_chain_mode` -- `"walls_only"` (default, proteins+rRNA) or `"assembly_all"` (everything in first assembly).
- `occupancy_exclude_trna` -- auto-detect and exclude tRNA chains. Default True.
- `occupancy_exclude_auth_asym_ids` -- tuple of chain IDs to manually exclude. Default empty.
- `cylinder_radius_A` -- radius of the search cylinder. Default 60.
- `cylinder_height_A` -- height of the search cylinder. Default 150.

### 40_empty_space
Voxelizes the search cylinder, marks voxels occupied by atoms, then extracts empty voxels. Empty voxels outside the exterior shell (Stage 20) are discarded.

Supports two occupancy backends per grid level:
- `legacy_kdtree` -- ball query around each voxel center
- `edt` -- Euclidean distance transform (faster for fine grids)

Config parameters:
- `grid_levels` -- list of `GridLevelConfig` entries. Each has:
  - `name` -- identifier, e.g. `"level_0"`
  - `voxel_size_A` -- voxel edge length. Default 1.0 for level_0.
  - `occupancy_backend` -- `"legacy_kdtree"` or `"edt"`. Default `"legacy_kdtree"`.
  - `uniform_atom_radius_A` -- radius around each atom center that counts as occupied. Default 2.0.

### 50_clustering
Two-pass DBSCAN on the empty voxels from Stage 40 to isolate the tunnel void from other cavities.

1. **Coarse pass**: large eps, high min_samples -- merges nearby regions
2. **Refine pass**: tighter eps on the coarse winner -- cleans edges

Cluster selection uses proximity to the constriction site (the cluster containing or nearest to the constriction point is the tunnel).

Optionally generates a coarse mesh (level_0) from the refined cluster.

Config parameters:
- `dbscan_level0_coarse_eps_A` -- coarse DBSCAN neighborhood radius. Default 5.5.
- `dbscan_level0_coarse_min_samples` -- coarse DBSCAN density threshold. Default 600.
- `dbscan_level0_refine_eps_A` -- refine DBSCAN neighborhood radius. Default 3.5.
- `dbscan_level0_refine_min_samples` -- refine DBSCAN density threshold. Default 175.
- `mesh_level0_enable` -- whether to generate a mesh at this stage. Default True.

### 55_grid_refine
Re-voxelizes a tight bounding box around the Stage 50 cluster at a finer resolution (typically 0.5A), recomputes occupancy via EDT, extracts boundary voxels, and runs another two-pass DBSCAN. Generates the high-detail mesh (level_1).

Config parameters:
- `refine_voxel_size_A` -- voxel size for the fine grid. Default 0.5.
- `refine_roi_pad_A` -- padding around Stage 50 cluster bbox. Default 10.
- `refine_atom_radius_A` -- atom occupancy radius on the fine grid. Default 2.0.
- `refine_keep_within_A` -- max distance from Stage 50 cluster to keep fine voxels. Default 6.0.
- `refine_occ_close_iters` -- morphological closing iterations on occupancy. Default 0.
- `refine_void_open_iters` -- morphological opening iterations on void mask (removes thin bridges). Default 1.
- `refine_forbid_roi_boundary` -- zero out voxels touching the ROI boundary. Default True.
- `dbscan_level1_coarse_eps_A` -- default 3.0.
- `dbscan_level1_coarse_min_samples` -- default 30.
- `dbscan_level1_refine_eps_A` -- default 3.0.
- `dbscan_level1_refine_min_samples` -- default 20.
- `refine_dbscan_max_points` -- cap on points sent to DBSCAN (0 = no cap). Default 0.
- `refine_dbscan_seed` -- RNG seed for subsampling. Default 0.
- `mesh_level1_enable` -- whether to generate the fine mesh. Default True.

### 70_mesh_validate
Picks the best available mesh (prefers level_1, falls back to level_0), validates watertightness, and copies final + comparison meshes to the run root.

No config parameters.

## Mesh outputs per run

All meshes are saved as both binary PLY and ASCII PLY (suffix `_ascii.ply`).
```
runs/{RCSB_ID}/{RUN_ID}/
  tunnel_mesh.ply                    # final mesh (binary)
  tunnel_mesh_ascii.ply              # final mesh (ASCII, for viewers)
  mesh_level_0.ply                   # coarse mesh (1.0A grid)
  mesh_level_0_ascii.ply
  mesh_level_0_pre_smooth.ply        # coarse mesh before Taubin smoothing
  mesh_level_0_pre_smooth_ascii.ply
  mesh_level_1.ply                   # fine mesh (0.5A grid)
  mesh_level_1_ascii.ply
  mesh_level_1_pre_smooth.ply        # fine mesh before Taubin smoothing
  mesh_level_1_pre_smooth_ascii.ply
  manifest.json                      # full run provenance
  stage/                             # per-stage artifacts
    ...
```

## Smoothing parameters

These control how the final mesh surface looks. If the mesh is too rough or too smooth, these are what to adjust.

`mesh_gaussian_sigma_voxels` (default 1.5) -- Gaussian blur applied to the binary volume *before* marching cubes. Higher values produce a smoother isosurface but lose fine detail. Range to try: 0.5 (sharp) to 3.0 (very smooth). This is the single most impactful smoothing parameter.

`mesh_level0_smooth_iters` (default 20) -- Taubin smoothing iterations on the coarse (level_0) mesh *after* marching cubes. More iterations = smoother. 0 disables post-MC smoothing entirely.

`mesh_level1_smooth_iters` (default 40) -- same, for the fine (level_1) mesh. Higher default because the 0.5A grid produces more surface detail that benefits from smoothing.

`mesh_smooth_method` (default `"taubin"`) -- smoothing algorithm. Taubin preserves volume better than Laplacian. Change to `"laplacian"` only if you want aggressive shrinkage.

`mesh_taubin_pass_band` (default 0.1) -- Taubin frequency cutoff. Lower = more aggressive smoothing per iteration. Range: 0.01 (aggressive) to 0.5 (gentle).

`mesh_fill_holes_A` (default 100.0) -- max hole size (in angstroms) to fill after marching cubes. Prevents small gaps from breaking watertightness.

`mesh_atom_clearance_A` (default 1.5) -- after smoothing, any mesh vertex closer than this to an atom is pushed outward. Prevents the smoothed surface from penetrating into occupied space. Set to 0 to disable.

### Quick recipes

**Smoother tunnel**: increase `mesh_gaussian_sigma_voxels` to 2.0-2.5 and/or increase `mesh_level1_smooth_iters` to 60-80.

**Sharper tunnel** (more surface detail): decrease `mesh_gaussian_sigma_voxels` to 0.8-1.0, decrease `mesh_level1_smooth_iters` to 10-20.

**Debug (no smoothing at all)**: set `mesh_gaussian_sigma_voxels` to 0, `mesh_level1_smooth_iters` to 0. The pre-smooth meshes are also saved automatically for comparison.

## Run naming

Runs are named `{SEQ}_{TIMESTAMP}_{HASH}`, e.g. `003_20260211_143022_fe15457ff09fb93b`. The sequential index increments per structure, so the latest run is always the highest number.