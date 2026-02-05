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
