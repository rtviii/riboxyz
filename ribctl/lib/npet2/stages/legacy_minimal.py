# ribctl/lib/npet2/stages/legacy_minimal.py
# npet2/stages/legacy_minimal.py

from __future__ import annotations
import json
from pathlib import Path
import time
from typing import Any, Dict, List, Tuple
import numpy as np
import pyvista as pv
import open3d as o3d

from npet2.backends.grid_occupancy import (
    connected_components_3d,
    occupancy_via_edt,
)
from npet2.backends.meshing import save_mesh_with_ascii
from npet2.core.cache import StageCacheKey
from npet2.core.pipeline import Stage
from npet2.core.ribosome_types import RibosomeProfile
from npet2.core.structure_selection import (
    intersect_with_first_assembly,
    ribosome_wall_auth_asym_ids,
    tunnel_debris_chains,
    atom_inclusion_policy,
)
from npet2.core.types import StageContext, ArtifactType

from scipy import ndimage

from npet2.backends.geometry import (
    cif_to_point_cloud,
    fast_normal_estimation,
    quick_surface_points,
    validate_mesh_pyvista,
    apply_poisson_reconstruction,
    filter_residues_parallel,
    transform_points_to_C0,
    transform_points_from_C0,
    create_point_cloud_mask,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    estimate_normals,
)
from npet2.stages.grid_refine import (
    _make_bbox_grid,
    _points_to_ijk,
    _valid_ijk,
    _voxel_centers_from_indices,
)

# ... rest of the file is identical, except within _generate_mesh methods,
# replace:
#   from ribctl.lib.npet.kdtree_approach import ...
#   from npet2.backends.meshing import ...
# with:
#   from npet2.backends.meshing import ...
# (transform functions already imported at top)


def _residues_from_chain_ids(structure, chain_ids: set[str]):
    model = structure[0]
    residues = []
    for cid in chain_ids:
        if cid not in model:
            continue
        chain = model[cid]
        for r in chain.get_residues():
            if len(getattr(r, "child_list", [])) == 0:
                continue
            residues.append(r)
    return residues


def _pick_tunnel_cluster(
    clusters: dict[int, list],
    constr: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Pick the cluster whose points are closest to the constriction site."""
    constr = np.asarray(constr, dtype=np.float32).reshape(1, 3)
    best_id = -1
    best_dist = float("inf")
    for cid, pts_list in clusters.items():
        if cid == -1:
            continue
        pts = np.asarray(pts_list, dtype=np.float32)
        if pts.shape[0] == 0:
            continue
        dists = np.linalg.norm(pts - constr, axis=1)
        min_dist = float(dists.min())
        if min_dist < best_dist:
            best_dist = min_dist
            best_id = cid
    if best_id == -1:
        raise ValueError("No valid clusters found")
    print(f"  [cluster_select] picked cluster {best_id} "
          f"(n={len(clusters[best_id]):,}, dist_to_constriction={best_dist:.1f}A)")
    return np.asarray(clusters[best_id], dtype=np.float32), best_id

def _get_biopython_structure(ctx: StageContext):
    """Get biopython structure from ctx, or parse mmcif on demand."""
    bs = ctx.inputs.get("biopython_structure")
    if bs is not None:
        return bs
    from Bio.PDB.MMCIFParser import FastMMCIFParser
    mmcif_path = ctx.require("mmcif_path")
    bs = FastMMCIFParser(QUIET=True).get_structure(ctx.rcsb_id, mmcif_path)
    ctx.inputs["biopython_structure"] = bs
    return bs


class Stage20ExteriorShell(Stage):
    key = "20_exterior_shell"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A"  : c.shell_voxel_size_A,
            "dilation_iters": c.shell_dilation_iters,
            "closing_iters" : c.shell_closing_iters,
            "gaussian_sigma": c.shell_gaussian_sigma,
            "smooth_iters"  : c.shell_smooth_iters,
            "fill_holes"    : c.shell_fill_holes,
        }

    def run(self, ctx: StageContext) -> None:
        stage_cache = ctx.require("stage_cache")
        inputs_fp = ctx.require("inputs_fp")
        params = self.params(ctx)

        key = StageCacheKey(
            stage=self.key,
            inputs_fp={"structure": inputs_fp["structure"]},
            params=params,
            impl_version="v2",
        )

        stage_dir = ctx.store.stage_dir(self.key)
        cached_files = [
            "alpha_shell.ply",
            "alpha_shell_quality.json",
            "ribosome_ptcloud.npy",
        ]

        if stage_cache.has(key, required=["alpha_shell.ply", "alpha_shell_quality.json"]):
            stage_cache.copy_into(key, stage_dir, cached_files)
            quality = json.loads((stage_dir / "alpha_shell_quality.json").read_text())
            ctx.inputs["alpha_shell_path"] = str(stage_dir / "alpha_shell.ply")
            ctx.inputs["alpha_shell_watertight"] = bool(quality.get("watertight", False))
            ctx.store.register_file(name="alpha_shell_mesh", stage=self.key,
                                    type=ArtifactType.PLY_MESH, path=stage_dir / "alpha_shell.ply")
            ctx.store.register_file(name="alpha_shell_quality", stage=self.key,
                                    type=ArtifactType.JSON, path=stage_dir / "alpha_shell_quality.json")
            return

        c = ctx.config
        profile: RibosomeProfile = ctx.require("profile")
        cifpath = Path(ctx.require("mmcif_path"))

        ptcloud_path = stage_dir / "ribosome_ptcloud.npy"
        mesh_path    = stage_dir / "alpha_shell.ply"
        quality_path = stage_dir / "alpha_shell_quality.json"

        wall = ribosome_wall_auth_asym_ids(
            profile,
            exclude_trna=bool(getattr(ctx.config, "occupancy_exclude_trna", True)),
            extra_exclude=tunnel_debris_chains(ctx.rcsb_id, profile),
        )
        wall = intersect_with_first_assembly(profile, wall)
        ptcloud = cif_to_point_cloud(str(cifpath), sorted(wall), do_atoms=True)

        np.save(ptcloud_path, ptcloud)
        ctx.store.register_file(name="ribosome_ptcloud", stage=self.key,
                                type=ArtifactType.NUMPY, path=ptcloud_path)

        mesh, watertight = _make_voxel_shell(
            ptcloud,
            voxel_size_A   = c.shell_voxel_size_A,
            dilation_iters = c.shell_dilation_iters,
            closing_iters  = c.shell_closing_iters,
            gaussian_sigma = c.shell_gaussian_sigma,
            smooth_iters   = c.shell_smooth_iters,
            fill_holes     = c.shell_fill_holes,
        )

        mesh.save(str(mesh_path))

        quality = {
            "watertight" : bool(watertight),
            "n_points"   : int(mesh.n_points),
            "n_faces"    : int(mesh.n_faces_strict),
            "open_edges" : int(mesh.n_open_edges),
            "is_manifold": bool(mesh.is_manifold),
            "bounds"     : [float(x) for x in mesh.bounds],
        }
        quality_path.write_text(json.dumps(quality, indent=2))

        ctx.store.register_file(name="alpha_shell_mesh", stage=self.key,
                                type=ArtifactType.PLY_MESH, path=mesh_path)
        ctx.store.register_file(name="alpha_shell_quality", stage=self.key,
                                type=ArtifactType.JSON, path=quality_path)
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
        profile: RibosomeProfile = ctx.require("profile")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        policy = atom_inclusion_policy(profile, c, ctx.rcsb_id)

        occ_chain_ids = policy["wall_chain_ids"]
        seed_chain_ids = set(occ_chain_ids)

        structure = _get_biopython_structure(ctx)

        residues_seed = _residues_from_chain_ids(structure, seed_chain_ids)
        residues_occ = _residues_from_chain_ids(structure, occ_chain_ids)

        residues_seed = filter_residues_parallel(
            residues=residues_seed,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
        )
        residues_occ = filter_residues_parallel(
            residues=residues_occ,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
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

        out_seed = stage_dir / "region_atom_xyz.npy"
        np.save(out_seed, seed_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_seed,
            meta={
                "n": int(seed_points.shape[0]),
                "note": "seed atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz"] = seed_points

        out_occ = stage_dir / "region_atom_xyz_occ.npy"
        np.save(out_occ, occ_points)
        ctx.store.register_file(
            name="region_atom_xyz_occ",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_occ,
            meta={
                "n": int(occ_points.shape[0]),
                "note": "occupancy atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz_occ"] = occ_points

        policy_record = {
            "occupancy_chain_mode": getattr(c, "occupancy_chain_mode", "walls_only"),
            "exclude_trna": bool(getattr(c, "occupancy_exclude_trna", True)),
            "wall_chain_ids": sorted(occ_chain_ids),
            "excluded_chains": {k: v for k, v in policy["reasons"].items()},
            "policy_summary": (
                "INCLUDED: ribosomal proteins + rRNAs (with modified residues). "
                "EXCLUDED: waters, ions, nonpolymer ligands, tRNAs, debris chains."
            ),
        }
        (stage_dir / "atom_selection_policy.json").write_text(
            json.dumps(policy_record, indent=2)
        )

        print(
            f"[{self.key}] seed_atoms={seed_points.shape[0]:,} "
            f"occ_atoms={occ_points.shape[0]:,} occ_chains={len(occ_chain_ids)}"
        )
        if policy["reasons"]:
            excluded_summary = ", ".join(
                f"{k}({v})" for k, v in sorted(policy["reasons"].items())
            )
            print(f"[{self.key}] excluded: {excluded_summary}")

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

        from npet2.backends.grid_occupancy import (
            make_cylinder_grid,
            cylinder_mask,
            occupancy_via_edt,
            empty_points_from_mask,
            save_grid_npy,
            get_occupied_voxel_centers,
        )

        c = ctx.config

        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        # Use filtered atoms for clustering seed reference
        region_xyz_filtered = np.asarray(
            ctx.require("region_atom_xyz"), dtype=np.float32
        )

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
                    z_min=z_min,
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
                    z_min=z_min,
                )

                occ = occupancy_via_edt(
                    region_c0,
                    grid,
                    atom_radius_A=float(gl.uniform_atom_radius_A),
                )

                cyl2d = cylinder_mask(grid, radius_A=float(c.cylinder_radius_A))
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
    DBSCAN clustering on level_0 (coarse grid, typically 1.0A).
    
    Two-pass strategy:
      1. Coarse DBSCAN: merge regions, bridge gaps
      2. Refine DBSCAN: tighten on largest cluster from pass 1
    
    Cluster selection uses axial proximity to the PTC-Constriction axis
    rather than raw point count, which prevents the inter-subunit space
    from being picked over the actual tunnel.
    
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
            "cluster_selection": "axial_proximity",
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

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

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

        # Pick tunnel cluster from coarse (axial proximity, not largest)


        largest, largest_id = _pick_tunnel_cluster(clusters_coarse, constr)
        if largest.shape[0] == 0:
            raise ValueError(f"[{self.key}] tunnel cluster is empty")

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest)
        ctx.store.register_file(
            name="largest_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_largest,
            meta={"cluster_id": int(largest_id), "n": int(largest.shape[0]),
                  "selection": "axial_proximity"},
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

        # Pick tunnel cluster from refine


        refined, refined_id = _pick_tunnel_cluster(clusters_refine, constr)
        if refined.shape[0] == 0:
            raise ValueError(f"[{self.key}] refined cluster is empty")

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0]),
                  "selection": "axial_proximity"},
        )

        # Output for downstream
        ctx.inputs["refined_cluster"] = refined
        ctx.inputs["largest_cluster"] = largest

        print(f"[{self.key}] winner: coarse={largest.shape[0]:,} -> refine={refined.shape[0]:,}")

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
            if lab == -1:
                continue
            arr = np.asarray(plist, dtype=np.float32)
            if arr.size > 0:
                np.save(pass_dir / f"cluster_id{lab}.npy", arr)

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        import time
        from ribctl.lib.npet.kdtree_approach import transform_points_to_C0, transform_points_from_C0
        from npet2.backends.meshing import (
            mesh_from_binary_volume,
            voxelize_points,
            clip_mesh_to_atom_clearance,
            save_mesh_with_ascii,
        )

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        print(f"[{self.key}] generating mesh for {level_name}...")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        pts_c0 = transform_points_to_C0(points, ptc, constr).astype(np.float32)
        voxel = 1.0

        t0 = time.perf_counter()
        mask, origin = voxelize_points(pts_c0, voxel_size=voxel, pad_voxels=2)

        try:
            surf_c0, pre_smooth_c0 = mesh_from_binary_volume(
                mask, origin, voxel,
                gaussian_sigma_voxels=c.mesh_level0_gaussian_sigma,
                smooth_method=c.mesh_smooth_method,
                smooth_iters=c.mesh_level0_smooth_iters,
                taubin_pass_band=c.mesh_taubin_pass_band,
                fill_holes_size=c.mesh_fill_holes_A,
            )
        except ValueError as e:
            print(f"[{self.key}] MC mesh failed for {level_name}: {e}")
            return

        # Transform both meshes to world coordinates
        def _to_world(mesh_c0: pv.PolyData) -> pv.PolyData:
            pts_w = transform_points_from_C0(
                np.asarray(mesh_c0.points, dtype=np.float32), ptc, constr
            ).astype(np.float32)
            m = mesh_c0.copy(deep=True)
            m.points = pts_w
            return m

        pre_smooth_w = _to_world(pre_smooth_c0)
        pre_smooth_path = stage_dir / f"mesh_{level_name}_pre_smooth.ply"
        save_mesh_with_ascii(pre_smooth_w, pre_smooth_path, tag=f"{level_name}-pre-smooth")

        surf_w = _to_world(surf_c0)

        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)
        surf_w = clip_mesh_to_atom_clearance(surf_w, region_xyz, min_clearance_A=c.mesh_atom_clearance_A)

        dt = time.perf_counter() - t0
        print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
              f"{surf_w.n_faces:,} faces, watertight={surf_w.is_manifold and surf_w.n_open_edges == 0}")

        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        save_mesh_with_ascii(surf_w, mesh_path, tag=level_name)

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name, "method": "marching_cubes_taubin"},
        )
        print(f"[{self.key}] mesh saved: {mesh_path}")

class Stage70MeshValidate(Stage):
    key = "70_mesh_validate"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        import json
        import shutil

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

        # Try level_1 first (higher detail), fall back to level_0
        chosen_src = None
        chosen_label = None

        l1_path = ctx.inputs.get("level_1_mesh_path")
        if l1_path and Path(l1_path).exists():
            m = pv.read(l1_path)
            if m.is_manifold and m.n_open_edges == 0 and m.n_points > 0:
                chosen_src = l1_path
                chosen_label = "level_1"
                print(f"[{self.key}] using level_1 mesh (0.5A grid)")

        if chosen_src is None:
            # Look for level_0
            l0_path = ctx.store.run_dir / "stage" / "50_clustering" / "mesh_level_0.ply"
            if l0_path.exists():
                m = pv.read(str(l0_path))
                if m.n_points > 0:
                    chosen_src = str(l0_path)
                    chosen_label = "level_0"
                    print(f"[{self.key}] falling back to level_0 mesh (1.0A grid)")

        if chosen_src is None:
            raise ValueError(f"[{self.key}] no valid mesh found from any stage")

        shutil.copy2(chosen_src, mesh_path)
        final = pv.read(str(mesh_path))

        # Save final mesh as both binary and ASCII
        save_mesh_with_ascii(final, mesh_path, tag="final")

        st = _mesh_stats(final)
        watertight = final.is_manifold and final.n_open_edges == 0
        print(f"[{self.key}] final mesh: {st}, watertight={watertight}")

        if not watertight:
            raise ValueError(f"[{self.key}] final mesh is not watertight")

        ctx.store.register_file(
            name="tunnel_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"watertight": True, "source": chosen_label},
        )
        ctx.inputs["tunnel_mesh_path"] = str(mesh_path)

        # Also copy final meshes to run root for convenience
        root_mesh = ctx.store.run_dir / "tunnel_mesh.ply"
        root_mesh_ascii = ctx.store.run_dir / "tunnel_mesh_ascii.ply"
        shutil.copy2(str(mesh_path), str(root_mesh))
        ascii_src = mesh_path.parent / f"{mesh_path.stem}_ascii.ply"
        if ascii_src.exists():
            shutil.copy2(str(ascii_src), str(root_mesh_ascii))

        self._copy_comparison_meshes(ctx, stage_dir)

    def _copy_comparison_meshes(self, ctx, stage_dir):
        import shutil

        run_root = ctx.store.run_dir

        for stage_name, level, voxel in [
            ("50_clustering", "level_0", 1.0),
            ("55_grid_refine", "level_1", 0.5),
        ]:
            src_dir = ctx.store.run_dir / "stage" / stage_name

            # Post-smooth mesh
            for suffix in [f"mesh_{level}.ply", f"mesh_{level}_ascii.ply"]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    # Also put in run root
                    shutil.copy2(src, run_root / suffix)

            # Pre-smooth mesh
            for suffix in [
                f"mesh_{level}_pre_smooth.ply",
                f"mesh_{level}_pre_smooth_ascii.ply",
            ]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    shutil.copy2(src, run_root / suffix)

            if (src_dir / f"mesh_{level}.ply").exists():
                print(
                    f"[{self.key}]   copied {level} mesh ({voxel}A grid) + pre-smooth"
                )
