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

        # Filtered residues (for clustering/DBSCAN seed)
        residues_filtered = ribosome_entities(
            rcsb_id=ctx.rcsb_id,
            cifpath=cifpath,
            level="R",
            skip_nascent_chain=skip,
        )

        residues_filtered = filter_residues_parallel(
            residues=residues_filtered,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
        )

        filtered_points = np.asarray(
            [atom.get_coord() for r in residues_filtered for atom in r.child_list],
            dtype=np.float32,
        )

        out = ctx.store.stage_dir(self.key) / "region_atom_xyz.npy"
        np.save(out, filtered_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out,
            meta={"n": int(filtered_points.shape[0]), "note": "filtered atoms for clustering"},
        )

        ctx.inputs["region_atom_xyz"] = filtered_points

        # ALL residues (for occupancy to prevent mesh interference)
        residues_all = ribosome_entities(
            rcsb_id=ctx.rcsb_id,
            cifpath=cifpath,
            level="R",
            skip_nascent_chain=[],  # don't skip anything
        )

        residues_all = filter_residues_parallel(
            residues=residues_all,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
        )

        all_points = np.asarray(
            [atom.get_coord() for r in residues_all for atom in r.child_list],
            dtype=np.float32,
        )

        out_all = ctx.store.stage_dir(self.key) / "region_atom_xyz_all.npy"
        np.save(out_all, all_points)
        ctx.store.register_file(
            name="region_atom_xyz_all",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_all,
            meta={"n": int(all_points.shape[0]), "note": "ALL atoms for occupancy (prevents mesh interference)"},
        )

        ctx.inputs["region_atom_xyz_all"] = all_points
        
        print(f"[{self.key}] filtered atoms: {filtered_points.shape[0]:,}, all atoms: {all_points.shape[0]:,}")

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
        
        # Use ALL atoms for occupancy (prevents mesh interference)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        
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