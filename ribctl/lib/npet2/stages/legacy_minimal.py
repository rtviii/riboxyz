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
