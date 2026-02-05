# __scripts/npet2_viz_enhanced.py
#!/usr/bin/env python3
"""
NPET2 multi-panel visualization (coarse -> refine -> grid-refine DBSCAN -> final).

Panels (1x3):
  (0) Stage40 empty points + Stage50 coarse DBSCAN
  (1) Stage50 refine DBSCAN + Stage55 ROI bbox + Stage50 winners
  (2) Stage55 refined-grid DBSCAN (coarse/refine toggle) + Stage55 final refined surface + Stage70 mesh

New Stage55 DBSCAN artifacts expected (from the Stage55GridRefine05 replacement):
  stage/55_grid_refine/dbscan_coarse/points.npy + labels.npy
  stage/55_grid_refine/dbscan_refine/points.npy + labels.npy
  stage/55_grid_refine/refined_surface_points_level_1.npy

Usage examples:
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --shell
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --dbscan_show_noise --downsample 10
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --stage55_dbscan coarse
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --stage55_dbscan refine
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --show_voxel_surface
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Optional

import numpy as np
import pyvista as pv


# ----------------------------
# Basic IO helpers
# ----------------------------


def find_latest_run_dir(runs_root: Path, rcsb: str) -> Path:
    rdir = runs_root / rcsb.upper()
    if not rdir.exists():
        raise FileNotFoundError(f"No runs for {rcsb} under {runs_root}")
    candidates = [p for p in rdir.iterdir() if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No run subdirs for {rcsb} under {rdir}")
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def load_npy_points(path: Path) -> np.ndarray:
    a = np.load(path)
    a = np.asarray(a, dtype=np.float32)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"{path} is not (N,3): got {a.shape}")
    return a


def _safe_read_mesh(path: Path) -> Optional[pv.DataSet]:
    try:
        return pv.read(str(path))
    except Exception as e:
        print(f"[viz] failed to read mesh {path}: {e}")
        return None


# ----------------------------
# Downsampling (avoid striping artifacts)
# ----------------------------


def random_downsample(pts: np.ndarray, downsample: int, seed: int = 0):
    """
    Avoid deterministic "striping" when pts are grid-ordered.
    Returns (pts_ds, idx_ds) where idx_ds indexes original pts.
    """
    if downsample <= 1:
        return pts, None
    n = int(pts.shape[0])
    k = max(1, n // int(downsample))
    rng = np.random.default_rng(int(seed))
    idx = rng.choice(n, size=k, replace=False)
    idx.sort()
    return pts[idx], idx


# ----------------------------
# ROI (bbox) drawing
# ----------------------------


def _bbox_corners(lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    return np.array(
        [
            [lo[0], lo[1], lo[2]],
            [hi[0], lo[1], lo[2]],
            [hi[0], hi[1], lo[2]],
            [lo[0], hi[1], lo[2]],
            [lo[0], lo[1], hi[2]],
            [hi[0], lo[1], hi[2]],
            [hi[0], hi[1], hi[2]],
            [lo[0], hi[1], hi[2]],
        ],
        dtype=np.float32,
    )


def add_roi_bbox_from_stage55(
    plotter: pv.Plotter,
    roi_json_path: Path,
    *,
    color="red",
    line_width=3,
    label="ROI bbox",
):
    """
    ROI bbox is stored in C0 frame; convert corners to world using saved ptc/constr.
    """
    try:
        roi = json.loads(roi_json_path.read_text())
    except Exception as e:
        print(f"[viz] failed to read ROI json {roi_json_path}: {e}")
        return

    lo = np.asarray(roi["lo"], dtype=np.float32)
    hi = np.asarray(roi["hi"], dtype=np.float32)
    corners_c0 = _bbox_corners(lo, hi)

    ptc = np.asarray(roi["transform"]["ptc"], dtype=np.float32)
    constr = np.asarray(roi["transform"]["constriction"], dtype=np.float32)

    # NOTE: keep import local so script works even if ribctl isn't in env for some reason
    from ribctl.lib.npet.kdtree_approach import transform_points_from_C0

    corners_w = transform_points_from_C0(corners_c0, ptc, constr).astype(np.float32)

    edges = [
        [0, 1],
        [1, 2],
        [2, 3],
        [3, 0],
        [4, 5],
        [5, 6],
        [6, 7],
        [7, 4],
        [0, 4],
        [1, 5],
        [2, 6],
        [3, 7],
    ]
    for i, (a, b) in enumerate(edges):
        line = pv.Line(corners_w[a], corners_w[b])
        plotter.add_mesh(
            line,
            color=color,
            line_width=line_width,
            opacity=1.0,
            label=label if i == 0 else None,
        )


# ----------------------------
# Point cloud helpers
# ----------------------------


def add_points(
    plotter: pv.Plotter,
    pts: np.ndarray,
    *,
    label: str,
    point_size: int,
    opacity: float,
    color: Optional[str] = None,
):
    if pts.size == 0:
        print(f"[viz] {label}: empty array")
        return
    poly = pv.PolyData(pts)
    kwargs = dict(
        render_points_as_spheres=True,
        point_size=int(point_size),
        opacity=float(opacity),
        label=label,
    )
    if color is not None:
        kwargs["color"] = color
    plotter.add_points(poly, **kwargs)


def add_dbscan_pass(
    plotter: pv.Plotter,
    pass_dir: Path,
    *,
    label: str,
    point_size: int,
    downsample: int,
    opacity: float,
    show_noise: bool,
    seed: int = 0,
):
    pts_p = pass_dir / "points.npy"
    lab_p = pass_dir / "labels.npy"
    if not (pts_p.exists() and lab_p.exists()):
        print(f"[viz] missing DBSCAN pass files in {pass_dir}")
        return

    pts = np.load(pts_p).astype(np.float32)
    lab = np.load(lab_p).astype(np.int32)

    if downsample > 1:
        pts, idx = random_downsample(pts, downsample, seed=seed)
        lab = lab[idx]

    # optionally draw noise separately
    if show_noise:
        m_noise = lab == -1
        if m_noise.any():
            add_points(
                plotter,
                pts[m_noise],
                label=f"{label} noise",
                point_size=max(1, point_size - 1),
                opacity=min(opacity, 0.25),
                color="gray",
            )

    # keep non-noise for colormap
    m = lab != -1
    pts, lab = pts[m], lab[m]
    if pts.shape[0] == 0:
        print(f"[viz] no non-noise points for {label}")
        return

    # shift labels to start at 0
    lab = lab - lab.min()

    poly = pv.PolyData(pts)
    poly["cluster"] = lab
    plotter.add_points(
        poly,
        scalars="cluster",
        cmap="tab20",
        render_points_as_spheres=True,
        point_size=int(point_size),
        opacity=float(opacity),
        label=label,
    )


# ----------------------------
# Voxel mask -> surface (optional)
# ----------------------------


def voxel_mask_to_surface(mask_path: Path, spec_path: Path) -> Optional[pv.PolyData]:
    """
    Deterministic marching-cubes style surface from a binary mask saved by Stage55.
    This returns a surface in WORLD coords.
    """
    if not (mask_path.exists() and spec_path.exists()):
        return None

    try:
        spec = json.loads(spec_path.read_text())
        vol = np.load(mask_path).astype(np.uint8)
    except Exception as e:
        print(f"[viz] failed to load voxel mask/spec: {e}")
        return None

    voxel = float(spec["voxel_size_A"])
    origin = np.asarray(spec["origin"], dtype=np.float32)
    ptc = np.asarray(spec["transform"]["ptc"], dtype=np.float32)
    constr = np.asarray(spec["transform"]["constriction"], dtype=np.float32)

    vol_pad = np.pad(vol, 1, constant_values=0)
    origin_pad = origin - voxel

    img = pv.ImageData(
        dimensions=(vol_pad.shape[0] + 1, vol_pad.shape[1] + 1, vol_pad.shape[2] + 1),
        spacing=(voxel, voxel, voxel),
        origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
    )
    img.cell_data["void"] = vol_pad.ravel(order="F")

    surf_c0 = img.contour(isosurfaces=[0.5], scalars="void").triangulate()
    if surf_c0.n_points == 0:
        return None

    from ribctl.lib.npet.kdtree_approach import transform_points_from_C0

    pts_w = transform_points_from_C0(
        np.asarray(surf_c0.points, dtype=np.float32), ptc, constr
    ).astype(np.float32)

    surf_w = surf_c0.copy(deep=True)
    surf_w.points = pts_w
    return surf_w


# ----------------------------
# Main
# ----------------------------


def main():
    ap = argparse.ArgumentParser(description="NPET2 multi-panel visualization")
    ap.add_argument("--runs_root", default="/Users/rtviii/dev/riboxyz/NPET2/runs")
    ap.add_argument("--rcsb", required=True)
    ap.add_argument("--run_dir", default=None)

    # toggles
    ap.add_argument(
        "--shell", action="store_true", help="Show alpha shell (faint) in all panels"
    )
    ap.add_argument(
        "--show_voxel_surface",
        action="store_true",
        help="Show Stage55 selected void component surface in panel 2",
    )

    # Stage55 DBSCAN selection (new)
    ap.add_argument(
        "--stage55_dbscan",
        choices=["off", "coarse", "refine"],
        default="refine",
        help="Overlay Stage55 refined-grid DBSCAN clusters in panel 2 (off/coarse/refine)",
    )

    # appearance
    ap.add_argument("--point_size", type=int, default=6)
    ap.add_argument(
        "--downsample",
        type=int,
        default=3,
        help="Random downsample factor (>=1). Higher = fewer points.",
    )
    ap.add_argument("--opacity", type=float, default=0.75)
    ap.add_argument("--dbscan_show_noise", action="store_true")
    ap.add_argument("--legend", action="store_true", help="Show legend (can clutter)")
    ap.add_argument("--seed", type=int, default=0, help="Seed for random downsampling")

    args = ap.parse_args()

    runs_root = Path(args.runs_root)
    run_dir = (
        Path(args.run_dir)
        if args.run_dir
        else find_latest_run_dir(runs_root, args.rcsb)
    )
    print(f"[viz] Using run: {run_dir}")

    stage = run_dir / "stage"
    st20 = stage / "20_exterior_shell"
    st40 = stage / "40_empty_space"
    st50 = stage / "50_clustering"
    st55 = stage / "55_grid_refine"
    st70 = stage / "70_mesh_validate"

    # Common assets
    alpha_shell = st20 / "alpha_shell.ply"
    roi_json = st55 / "roi_bbox_c0.json"
    mesh_path = st70 / "npet2_tunnel_mesh.ply"

    # DBSCAN pass dirs (Stage50)
    coarse_pass = st50 / "coarse"
    refine_pass = st50 / "refine"

    # Stage40 points (last grid level empty points; plus maybe level_0 explicitly)
    empty_l0 = st40 / "empty_points_level_0.npy"
    empty_last = (
        st40 / "empty_points_level_1.npy"
    )  # if you have 2 levels; ok if missing

    # Stage55 refined surface points (final selection)
    refined_surface = st55 / "refined_surface_points_level_1.npy"

    # Stage55 DBSCAN dirs (new Stage55 behavior)
    st55_db_coarse = st55 / "dbscan_coarse"
    st55_db_refine = st55 / "dbscan_refine"

    # Stage55 voxel mask surface (optional legacy/debug)
    sel_mask = st55 / "selected_void_component_mask_level_1.npy"
    grid_spec = st55 / "grid_spec_level_1.json"

    # ----------------------------
    # Plotter: 1 row x 3 columns
    # ----------------------------
    pl = pv.Plotter(shape=(1, 3), border=True, window_size=(2100, 750))
    pl.set_background("white")

    def _maybe_add_shell(p: pv.Plotter):
        if not args.shell:
            return
        if not alpha_shell.exists():
            print(f"[viz] missing alpha shell: {alpha_shell}")
            return
        m = _safe_read_mesh(alpha_shell)
        if m is None:
            return
        p.add_mesh(m, opacity=0.08, color="lightgray", label="alpha_shell")

    def _maybe_add_roi(p: pv.Plotter):
        if not roi_json.exists():
            print(f"[viz] missing ROI json: {roi_json}")
            return
        add_roi_bbox_from_stage55(
            p, roi_json, color="red", line_width=3, label="ROI bbox"
        )

    # -------- Panel 0: Stage40 + Stage50 coarse
    pl.subplot(0, 0)
    pl.add_text("Panel 0: Stage40 + DBSCAN coarse", font_size=12)
    _maybe_add_shell(pl)
    _maybe_add_roi(pl)

    if empty_l0.exists():
        pts = load_npy_points(empty_l0)
        pts, _ = random_downsample(pts, args.downsample, seed=args.seed)
        add_points(
            pl,
            pts,
            label="empty_points L0",
            point_size=args.point_size,
            opacity=min(args.opacity, 0.5),
            color="steelblue",
        )
    else:
        if empty_last.exists():
            pts = load_npy_points(empty_last)
            pts, _ = random_downsample(pts, args.downsample, seed=args.seed)
            add_points(
                pl,
                pts,
                label="empty_points",
                point_size=args.point_size,
                opacity=min(args.opacity, 0.5),
                color="steelblue",
            )
        else:
            print(f"[viz] missing empty points: {empty_l0} and {empty_last}")

    if coarse_pass.exists():
        add_dbscan_pass(
            pl,
            coarse_pass,
            label="DBSCAN coarse",
            point_size=args.point_size,
            downsample=args.downsample,
            opacity=args.opacity,
            show_noise=args.dbscan_show_noise,
            seed=args.seed,
        )
    else:
        print(f"[viz] missing coarse pass dir: {coarse_pass}")

    pl.add_axes()
    pl.show_grid()

    # -------- Panel 1: Stage50 refine + ROI + winners
    pl.subplot(0, 1)
    pl.add_text("Panel 1: DBSCAN refine + ROI", font_size=12)
    _maybe_add_shell(pl)
    _maybe_add_roi(pl)

    if refine_pass.exists():
        add_dbscan_pass(
            pl,
            refine_pass,
            label="DBSCAN refine",
            point_size=args.point_size,
            downsample=args.downsample,
            opacity=args.opacity,
            show_noise=args.dbscan_show_noise,
            seed=args.seed,
        )
    else:
        print(f"[viz] missing refine pass dir: {refine_pass}")

    # Winners for reference
    p_largest = st50 / "largest_cluster.npy"
    p_refined = st50 / "refined_cluster.npy"

    if p_largest.exists():
        pts = load_npy_points(p_largest)
        pts, _ = random_downsample(pts, args.downsample, seed=args.seed)
        add_points(
            pl,
            pts,
            label="largest_cluster",
            point_size=args.point_size + 1,
            opacity=1.0,
            color="black",
        )
    else:
        print(f"[viz] missing {p_largest}")

    if p_refined.exists():
        pts = load_npy_points(p_refined)
        pts, _ = random_downsample(pts, args.downsample, seed=args.seed)
        add_points(
            pl,
            pts,
            label="refined_cluster",
            point_size=args.point_size + 2,
            opacity=1.0,
            color="darkred",
        )
    else:
        print(f"[viz] missing {p_refined}")

    pl.add_axes()
    pl.show_grid()

    # -------- Panel 2: Stage55 DBSCAN + Stage55 final surface + Stage70 mesh
    pl.subplot(0, 2)
    pl.add_text("Panel 2: Stage55 DBSCAN + surface + mesh", font_size=12)
    _maybe_add_shell(pl)
    _maybe_add_roi(pl)

    # Stage55 DBSCAN overlay (new)
    if args.stage55_dbscan != "off":
        pass_dir = st55_db_refine if args.stage55_dbscan == "refine" else st55_db_coarse
        if pass_dir.exists():
            add_dbscan_pass(
                pl,
                pass_dir,
                label=f"Stage55 DBSCAN {args.stage55_dbscan}",
                point_size=args.point_size,
                downsample=args.downsample,
                opacity=args.opacity,
                show_noise=args.dbscan_show_noise,
                seed=args.seed,
            )
        else:
            print(f"[viz] missing Stage55 DBSCAN dir: {pass_dir}")

    # Stage55 final selected surface points
    if refined_surface.exists():
        pts = load_npy_points(refined_surface)
        pts, _ = random_downsample(pts, args.downsample, seed=args.seed)
        add_points(
            pl,
            pts,
            label="stage55 refined_surface_points",
            point_size=args.point_size,
            opacity=min(1.0, max(0.35, args.opacity)),
            color="royalblue",
        )
    else:
        print(f"[viz] missing refined surface points: {refined_surface}")

    # Optional voxel surface (legacy/debug; may be absent with DBSCAN-based Stage55)
    if args.show_voxel_surface:
        surf = voxel_mask_to_surface(sel_mask, grid_spec)
        if surf is None:
            print(f"[viz] missing voxel surface inputs: {sel_mask} / {grid_spec}")
        else:
            pl.add_mesh(
                surf,
                opacity=0.20,
                color="orange",
                label="void component surface (voxel)",
            )

    # Final mesh
    if mesh_path.exists():
        m = _safe_read_mesh(mesh_path)
        if m is not None:
            pl.add_mesh(m, opacity=0.35, color="tan", label="final tunnel mesh")
    else:
        print(f"[viz] missing final mesh: {mesh_path}")

    pl.add_axes()
    pl.show_grid()

    # shared camera (so all panels comparable)
    try:
        pl.link_views()
    except Exception:
        pass

    if args.legend:
        pl.add_legend(bcolor="white", size=(0.25, 0.25))

    pl.show()


if __name__ == "__main__":
    main()
