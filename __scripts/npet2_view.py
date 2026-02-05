#!/usr/bin/env python3
from __future__ import annotations

import argparse
import glob
import json
import os
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional

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


def parse_cluster_id_from_name(path: Path) -> int:
    """
    Extract ..._id<cid>.npy
    """
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
):
    cloud = pv.PolyData(pts)
    plotter.add_points(
        cloud,
        render_points_as_spheres=True,
        point_size=point_size,
        opacity=opacity,
        label=label,
    )


def add_roi_bbox(plotter: pv.Plotter, roi_json_path: Path, *, label: str = "ROI bbox"):
    roi = json.loads(roi_json_path.read_text())
    lo = np.array(roi["lo"], dtype=np.float32)
    hi = np.array(roi["hi"], dtype=np.float32)
    bounds = (lo[0], hi[0], lo[1], hi[1], lo[2], hi[2])
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(box, style="wireframe", line_width=2, color="black", label=label)


def add_mesh(
    plotter: pv.Plotter, mesh_path: Path, label: str, *, opacity: float = 0.15
):
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
    """
    Combine all cluster files into one PolyData with per-point scalar cluster_id.
    This guarantees distinct colors even if many files exist.
    """
    if not paths:
        return

    all_pts: List[np.ndarray] = []
    all_lbl: List[np.ndarray] = []

    # Map cluster ids to a dense 0..K-1 palette index for prettier categorical coloring
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

    # tab20 is a good categorical map for <=20 clusters; if you have more, it still cycles.
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
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs_root", default="/Users/rtviii/dev/riboxyz/NPET2/runs")
    ap.add_argument("--rcsb", required=True)
    ap.add_argument(
        "--run_dir",
        default=None,
        help="Optional explicit run dir; otherwise picks latest for RCSB",
    )

    ap.add_argument("--shell", action="store_true")
    ap.add_argument("--tunnel_mesh", action="store_true")
    ap.add_argument("--roi", action="store_true")

    ap.add_argument("--points0", action="store_true")
    ap.add_argument("--points1", action="store_true")

    ap.add_argument(
        "--stage40_coarse",
        action="store_true",
        help="Show Stage40 coarse_cluster_*.npy",
    )
    ap.add_argument(
        "--stage50_coarse",
        action="store_true",
        help="Show Stage50 coarse_cluster_*.npy",
    )
    ap.add_argument(
        "--stage50_refine",
        action="store_true",
        help="Show Stage50 refine_cluster_*.npy",
    )

    ap.add_argument("--largest", action="store_true")
    ap.add_argument("--refined", action="store_true")

    ap.add_argument("--point_size", type=int, default=3)
    ap.add_argument("--downsample", type=int, default=1)
    ap.add_argument("--opacity", type=float, default=1.0)

    args = ap.parse_args()

    runs_root = Path(args.runs_root)
    if args.run_dir:
        run_dir = Path(args.run_dir)
    else:
        run_dir = find_latest_run_dir(runs_root, args.rcsb)

    stage = run_dir / "stage"
    st20 = stage / "20_exterior_shell"
    st40 = stage / "40_empty_space"
    st50 = stage / "50_clustering"
    st70 = stage / "70_mesh_validate"

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
            add_roi_bbox(pl, roi_path, label="ROI bbox")

    # Empty points
    if args.points0:
        p = st40 / "empty_points_level_0.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[:: args.downsample]
            add_points(
                pl,
                pts,
                "empty_points_level_0",
                point_size=args.point_size,
                opacity=args.opacity,
            )

    if args.points1:
        p = st40 / "empty_points_level_1.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[:: args.downsample]
            add_points(
                pl,
                pts,
                "empty_points_level_1",
                point_size=args.point_size,
                opacity=args.opacity,
            )

    # Cluster artifacts (colored by cluster scalar)
    if args.stage40_coarse:
        paths = sorted(
            Path(p).resolve() for p in glob.glob(str(st40 / "coarse_cluster_*_id*.npy"))
        )
        add_cluster_group(
            pl,
            paths,
            "Stage40 coarse clusters",
            point_size=args.point_size,
            downsample=args.downsample,
            opacity=args.opacity,
        )

    if args.stage50_coarse:
        paths = sorted(
            Path(p).resolve() for p in glob.glob(str(st50 / "coarse_cluster_*_id*.npy"))
        )
        add_cluster_group(
            pl,
            paths,
            "Stage50 coarse clusters",
            point_size=args.point_size,
            downsample=args.downsample,
            opacity=args.opacity,
        )

    if args.stage50_refine:
        paths = sorted(
            Path(p).resolve() for p in glob.glob(str(st50 / "refine_cluster_*_id*.npy"))
        )
        add_cluster_group(
            pl,
            paths,
            "Stage50 refine clusters",
            point_size=args.point_size,
            downsample=args.downsample,
            opacity=args.opacity,
        )

    # Winners (single clouds; these can share colors if both enabled — that’s fine)
    if args.largest:
        p = st50 / "largest_cluster.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[:: args.downsample]
            add_points(
                pl, pts, "largest_cluster", point_size=args.point_size + 1, opacity=1.0
            )

    if args.refined:
        p = st50 / "refined_cluster.npy"
        if p.exists():
            pts = load_npy(p)
            if args.downsample > 1:
                pts = pts[:: args.downsample]
            add_points(
                pl, pts, "refined_cluster", point_size=args.point_size + 2, opacity=1.0
            )

    pl.add_axes()
    pl.show_grid()
    pl.add_legend(bcolor="white")
    pl.show()


if __name__ == "__main__":
    main()
