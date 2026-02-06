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


class Stage55GridRefine(Stage):
    key = "55_grid_refine"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A": float(c.refine_voxel_size_A),
            "roi_pad_A": float(c.refine_roi_pad_A),
            "atom_radius_A": float(c.refine_atom_radius_A),
            "keep_within_A": float(c.refine_keep_within_A),
            "occ_close_iters": int(c.refine_occ_close_iters),
            "void_open_iters": int(c.refine_void_open_iters),
            "forbid_roi_boundary": bool(c.refine_forbid_roi_boundary),
            "coarse_eps_A": float(c.dbscan_level1_coarse_eps_A),
            "coarse_min_samples": int(c.dbscan_level1_coarse_min_samples),
            "refine_eps_A": float(c.dbscan_level1_refine_eps_A),
            "refine_min_samples": int(c.dbscan_level1_refine_min_samples),
            "dbscan_max_points": int(c.refine_dbscan_max_points),
            "dbscan_seed": int(c.refine_dbscan_seed),
            "mesh_enable": bool(getattr(c, "mesh_level1_enable", True)),
            "mesh_poisson_depth": int(getattr(c, "mesh_level1_poisson_depth", 8)),
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
        
        # Use ALL atoms for occupancy (prevents mesh interference)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        
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
        print(f"[{self.key}]   normal estimation: {dt1:.2f}s, {len(pcd.points):,} points")

        normals_path = stage_dir / f"normals_{level_name}.ply"
        try:
            o3d.io.write_point_cloud(str(normals_path), pcd)
        except Exception as e:
            print(f"[{self.key}] failed to write normals PCD for {level_name}: {e}")
            return

        mesh_path = stage_dir / f"mesh_{level_name}.ply"

        # Open3D's built-in Poisson (won't crash-terminate on failure)
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

        # Trim low-density vertices (Poisson extrapolates far from input points)
        densities = np.asarray(densities)
        density_threshold = np.quantile(densities, 0.01)
        vertices_to_remove = densities < density_threshold
        mesh_o3d.remove_vertices_by_mask(vertices_to_remove)

        o3d.io.write_triangle_mesh(str(mesh_path), mesh_o3d)

        # Cleanup in pyvista
        try:
            mesh = pv.read(str(mesh_path))
            mesh = mesh.fill_holes(2000.0)
            mesh = mesh.triangulate()
            mesh = mesh.connectivity(largest=True)  # remove floating pockets last
            mesh.save(str(mesh_path))

            # ASCII for pymol
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