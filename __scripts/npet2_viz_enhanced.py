# __scripts/npet2_viz_enhanced.py
#!/usr/bin/env python3
"""
NPET2 "pipeline at a glance" viewer (PyVista).

What it shows (side-by-side subplots in one window):
  [1] Stage50 coarse DBSCAN clusters (colored), optional noise
  [2] Stage50 refine DBSCAN clusters (colored), optional noise
  [3] Final: refine clusters + final tunnel mesh overlay (optional)

Also (in each subplot, when available):
  - alpha shell mesh (context)
  - ROI wireframe box (context) if roi_bbox*.json exists

Cluster labeling improvements:
  - annotates the TOP-N clusters (by size from index.json) with text at cluster centroids
  - prints a top-cluster summary box per subplot (ids + sizes)
  - uses discrete RGB colors (stable within a subplot) rather than a continuous scalar bar

If any artifact is missing, it prints a note and skips it.

Expected NPET2 Stage50 layout (as produced by the updated Stage50Clustering):
  stage/50_clustering/coarse/{points.npy, labels.npy, index.json, cluster_id<cid>.npy...}
  stage/50_clustering/refine/{points.npy, labels.npy, index.json, cluster_id<cid>.npy...}

Run:
  p3 __scripts/npet2_viz_enhanced.py --rcsb 7K00 --show_shell --show_mesh --label_top 8 --downsample_coarse 2 --downsample_refine 1

Tip:
  If rendering is slow, increase downsample_* or reduce point sizes / turn off render_points_as_spheres.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pyvista as pv


# ---------------------------
# Utilities
# ---------------------------


def find_latest_run_dir(runs_root: Path, rcsb: str) -> Path:
    rdir = runs_root / rcsb.upper()
    if not rdir.exists():
        raise FileNotFoundError(f"No runs for {rcsb} under {runs_root}")
    candidates = [p for p in rdir.iterdir() if p.is_dir()]
    if not candidates:
        raise FileNotFoundError(f"No run subdirs for {rcsb} under {rdir}")
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def read_json(path: Path) -> dict:
    return json.loads(path.read_text())


def safe_exists(path: Path, what: str) -> bool:
    if not path.exists():
        print(f"[viz] missing {what}: {path}")
        return False
    return True


def load_points_npy(path: Path) -> np.ndarray:
    a = np.load(path)
    a = np.asarray(a, dtype=np.float32)
    if a.ndim != 2 or a.shape[1] != 3:
        raise ValueError(f"{path} is not (N,3): got {a.shape}")
    return a


def load_landmarks(stage10_dir: Path) -> tuple[np.ndarray, np.ndarray]:
    ptc_json = stage10_dir / "ptc.json"
    constr_json = stage10_dir / "constriction_site.json"
    if not ptc_json.exists() or not constr_json.exists():
        raise FileNotFoundError(f"Landmarks not found under {stage10_dir}")
    ptc = np.array(read_json(ptc_json)["location"], dtype=np.float32)
    constr = np.array(read_json(constr_json)["location"], dtype=np.float32)
    return ptc, constr


def find_roi_json(stage40_dir: Path) -> Optional[Path]:
    """
    Best-effort: find any ROI bbox JSON. We accept:
      - roi_bbox_c0.json
      - roi_bbox*.json
    """
    p = stage40_dir / "roi_bbox_c0.json"
    if p.exists():
        return p
    cand = sorted(stage40_dir.glob("roi_bbox*.json"))
    return cand[0] if cand else None


def find_alpha_shell(stage20_dir: Path) -> Optional[Path]:
    p = stage20_dir / "alpha_shell.ply"
    return p if p.exists() else None


def find_final_mesh(stage70_dir: Path) -> Optional[Path]:
    # current Stage70 filename
    p = stage70_dir / "npet2_tunnel_mesh.ply"
    if p.exists():
        return p
    # fallback names (if you rename later)
    for alt in ("tunnel_mesh.ply", "final_mesh.ply", "npet_tunnel_mesh.ply"):
        q = stage70_dir / alt
        if q.exists():
            return q
    return None


# ---------------------------
# Drawing helpers
# ---------------------------


def add_mesh(
    plotter: pv.Plotter, mesh_path: Path, *, opacity: float, label: str
) -> None:
    try:
        mesh = pv.read(str(mesh_path))
        plotter.add_mesh(mesh, opacity=opacity, label=label, smooth_shading=True)
    except Exception as e:
        print(f"[viz] failed to read mesh {mesh_path}: {e}")


def add_roi_bbox(
    plotter: pv.Plotter,
    roi_json_path: Path,
    *,
    fallback_ptc: np.ndarray,
    fallback_constr: np.ndarray,
    color: str = "red",
    line_width: int = 4,
    label: str = "ROI",
) -> None:
    """
    Draw ROI bounding box. Accepts the historical structure:
      {
        "lo": [x,y,z], "hi": [x,y,z],
        "transform": {"ptc": [...], "constriction": [...]}
      }
    If transform is missing, uses fallback ptc/constr.
    Assumes lo/hi are in C0 coords; transforms to world.
    """
    try:
        from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
    except Exception as e:
        print(f"[viz] cannot import transform_points_from_C0: {e}")
        return

    try:
        roi = read_json(roi_json_path)
        lo = np.array(roi["lo"], dtype=np.float32)
        hi = np.array(roi["hi"], dtype=np.float32)

        if (
            isinstance(roi.get("transform"), dict)
            and "ptc" in roi["transform"]
            and "constriction" in roi["transform"]
        ):
            ptc = np.array(roi["transform"]["ptc"], dtype=np.float32)
            constr = np.array(roi["transform"]["constriction"], dtype=np.float32)
        else:
            ptc, constr = fallback_ptc, fallback_constr

        corners_c0 = np.array(
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

        corners_w = transform_points_from_C0(corners_c0, ptc, constr)

        edges = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 0),
            (4, 5),
            (5, 6),
            (6, 7),
            (7, 4),
            (0, 4),
            (1, 5),
            (2, 6),
            (3, 7),
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
    except Exception as e:
        print(f"[viz] failed ROI bbox {roi_json_path}: {e}")


def _compute_cluster_color_table(n: int) -> np.ndarray:
    """
    Return (n,3) uint8 RGB colors. Uses a matplotlib categorical colormap if available,
    else falls back to a simple hash-based coloring.
    """
    try:
        import matplotlib.cm as cm

        cmap = cm.get_cmap("tab20", max(n, 1))
        cols = (np.array([cmap(i)[:3] for i in range(n)]) * 255.0).astype(np.uint8)
        return cols
    except Exception:
        # deterministic-ish fallback
        cols = np.zeros((n, 3), dtype=np.uint8)
        for i in range(n):
            x = (i + 1) * 2654435761  # Knuth multiplicative hash
            cols[i] = np.array(
                [(x >> 0) & 255, (x >> 8) & 255, (x >> 16) & 255], dtype=np.uint8
            )
        return cols


def build_dbscan_rgb_polydata(
    points: np.ndarray,
    labels: np.ndarray,
    *,
    downsample: int,
    show_noise: bool,
    noise_downsample: int,
) -> tuple[Optional[pv.PolyData], Optional[pv.PolyData], Dict[int, int]]:
    """
    Build RGB-colored polydata for non-noise clusters.
    Returns:
      (clusters_poly, noise_poly, cid_to_color_index)

    cid_to_color_index maps original DBSCAN cluster_id -> color table index.
    """
    if points.shape[0] == 0:
        return None, None, {}

    pts = points
    lab = labels.astype(np.int32, copy=False)

    # Identify cluster ids excluding noise
    cluster_ids = sorted(int(x) for x in np.unique(lab) if int(x) != -1)
    cid_to_idx = {cid: i for i, cid in enumerate(cluster_ids)}
    colors = _compute_cluster_color_table(len(cluster_ids))

    # Downsample for plotting
    if downsample > 1:
        pts_ds = pts[::downsample]
        lab_ds = lab[::downsample]
    else:
        pts_ds = pts
        lab_ds = lab

    # Build cluster polydata
    m = lab_ds != -1
    pts_c = pts_ds[m]
    lab_c = lab_ds[m]

    clusters_poly = None
    if pts_c.shape[0] > 0:
        rgb = np.zeros((pts_c.shape[0], 3), dtype=np.uint8)
        # vectorized map: cluster id -> color index
        # (lab_c may be non-contiguous; use cid_to_idx)
        for cid, idx in cid_to_idx.items():
            mask = lab_c == cid
            if mask.any():
                rgb[mask] = colors[idx]
        clusters_poly = pv.PolyData(pts_c)
        clusters_poly["rgb"] = rgb

    # Build noise polydata
    noise_poly = None
    if show_noise:
        m_noise = lab_ds == -1
        noise_pts = pts_ds[m_noise]
        if noise_downsample > 1 and noise_pts.shape[0] > 0:
            noise_pts = noise_pts[::noise_downsample]
        if noise_pts.shape[0] > 0:
            noise_poly = pv.PolyData(noise_pts)

    return clusters_poly, noise_poly, cid_to_idx


def top_clusters_from_index(index: dict, top_n: int) -> list[dict]:
    clusters = [c for c in index.get("clusters", []) if int(c.get("id", -1)) != -1]
    clusters.sort(key=lambda c: int(c.get("n", 0)), reverse=True)
    return clusters[: max(0, int(top_n))]


def add_cluster_labels(
    plotter: pv.Plotter,
    pass_dir: Path,
    *,
    prefix: str,
    top_n: int,
    font_size: int,
    text_color: str,
) -> None:
    """
    Annotate TOP-N clusters with a text label at centroid.
    Uses index.json (if present) and cluster_id<cid>.npy files (fast + accurate).
    """
    idx_path = pass_dir / "index.json"
    if not idx_path.exists():
        print(f"[viz] missing index.json for labels: {idx_path}")
        return

    try:
        index = read_json(idx_path)
        top = top_clusters_from_index(index, top_n)
        if not top:
            return

        pts_for_labels = []
        txt_for_labels = []

        for c in top:
            cid = int(c["id"])
            n = int(c.get("n", 0))
            cpath = pass_dir / c.get("path", f"cluster_id{cid}.npy")
            if not cpath.exists():
                print(f"[viz] missing cluster points for label: {cpath}")
                continue
            pts = load_points_npy(cpath)
            if pts.shape[0] == 0:
                continue
            centroid = pts.mean(axis=0)
            pts_for_labels.append(centroid)
            txt_for_labels.append(f"{prefix}:{cid} (n={n})")

        if pts_for_labels:
            plotter.add_point_labels(
                np.array(pts_for_labels, dtype=np.float32),
                txt_for_labels,
                font_size=font_size,
                text_color=text_color,
                point_size=0,  # don't draw points for labels
                shape=None,
                show_points=False,
                always_visible=True,
            )
    except Exception as e:
        print(f"[viz] failed to add cluster labels for {pass_dir}: {e}")


def add_top_summary_text(
    plotter: pv.Plotter,
    pass_dir: Path,
    *,
    title: str,
    top_n: int,
    font_size: int,
) -> None:
    """
    Adds a small text box showing top cluster ids and sizes.
    """
    idx_path = pass_dir / "index.json"
    if not idx_path.exists():
        plotter.add_text(
            f"{title}\n(no index.json)", position="upper_right", font_size=font_size
        )
        return

    try:
        index = read_json(idx_path)
        top = top_clusters_from_index(index, top_n)
        lines = [title]
        if not top:
            lines.append("(no non-noise clusters)")
        else:
            for c in top:
                lines.append(f"id {int(c['id'])}: n={int(c['n'])}")
        plotter.add_text("\n".join(lines), position="upper_right", font_size=font_size)
    except Exception as e:
        plotter.add_text(
            f"{title}\n(index read failed)", position="upper_right", font_size=font_size
        )
        print(f"[viz] failed reading {idx_path}: {e}")


def add_dbscan_pass_subplot(
    plotter: pv.Plotter,
    pass_dir: Path,
    *,
    title: str,
    point_size: int,
    opacity: float,
    downsample: int,
    show_noise: bool,
    noise_point_size: int,
    noise_opacity: float,
    noise_downsample: int,
    render_spheres: bool,
    label_top: int,
    label_font_size: int,
    label_text_color: str,
    summary_top: int,
    summary_font_size: int,
) -> None:
    """
    Draw one DBSCAN pass (coarse/refine) into the current subplot.
    """
    pts_p = pass_dir / "points.npy"
    lab_p = pass_dir / "labels.npy"
    if not (pts_p.exists() and lab_p.exists()):
        print(
            f"[viz] missing DBSCAN pass files in {pass_dir} (need points.npy + labels.npy)"
        )
        plotter.add_text(
            f"{title}\n(missing artifacts)", position="upper_left", font_size=12
        )
        return

    pts = np.load(pts_p).astype(np.float32)
    lab = np.load(lab_p).astype(np.int32)

    clusters_poly, noise_poly, _ = build_dbscan_rgb_polydata(
        pts,
        lab,
        downsample=downsample,
        show_noise=show_noise,
        noise_downsample=noise_downsample,
    )

    plotter.add_text(title, position="upper_left", font_size=14)

    if clusters_poly is not None and clusters_poly.n_points > 0:
        plotter.add_points(
            clusters_poly,
            rgb=True,
            scalars="rgb",
            render_points_as_spheres=render_spheres,
            point_size=point_size,
            opacity=opacity,
            label=title,
        )
    else:
        print(
            f"[viz] {pass_dir}: no non-noise points to render (after downsampling/filtering)"
        )

    if noise_poly is not None and noise_poly.n_points > 0:
        plotter.add_points(
            noise_poly,
            color="gray",
            render_points_as_spheres=render_spheres,
            point_size=noise_point_size,
            opacity=noise_opacity,
            label=f"{title} noise",
        )

    # Labels (top clusters at centroid)
    add_cluster_labels(
        plotter,
        pass_dir,
        prefix=title.split()[-1].upper(),  # "COARSE"/"REFINE" style
        top_n=label_top,
        font_size=label_font_size,
        text_color=label_text_color,
    )

    # Summary text box
    add_top_summary_text(
        plotter,
        pass_dir,
        title="Top clusters",
        top_n=summary_top,
        font_size=summary_font_size,
    )


# ---------------------------
# Main
# ---------------------------


def main() -> None:
    ap = argparse.ArgumentParser(
        description="NPET2 pipeline-at-a-glance viewer (PyVista)"
    )

    ap.add_argument("--runs_root", default="/Users/rtviii/dev/riboxyz/NPET2/runs")
    ap.add_argument("--rcsb", required=True)
    ap.add_argument(
        "--run_dir",
        default=None,
        help="If omitted, uses latest run under runs_root/RCsb",
    )

    # Context overlays
    ap.add_argument(
        "--show_shell",
        action="store_true",
        help="Show alpha shell in each subplot (if present)",
    )
    ap.add_argument("--shell_opacity", type=float, default=0.10)

    ap.add_argument(
        "--show_roi",
        action="store_true",
        help="Show ROI wireframe in each subplot (if roi_bbox*.json exists)",
    )
    ap.add_argument("--roi_color", default="red")
    ap.add_argument("--roi_line_width", type=int, default=5)

    ap.add_argument(
        "--show_mesh",
        action="store_true",
        help="Show final tunnel mesh in the LAST subplot (if present)",
    )
    ap.add_argument("--mesh_opacity", type=float, default=0.35)

    # Appearance / performance
    ap.add_argument(
        "--render_spheres",
        action="store_true",
        help="Render points as spheres (prettier, slower)",
    )
    ap.add_argument(
        "--link_views", action="store_true", help="Link cameras across subplots"
    )

    # DBSCAN rendering params (coarse / refine)
    ap.add_argument("--downsample_coarse", type=int, default=2)
    ap.add_argument("--downsample_refine", type=int, default=1)

    ap.add_argument("--point_size_coarse", type=int, default=6)
    ap.add_argument("--point_size_refine", type=int, default=7)
    ap.add_argument("--opacity_coarse", type=float, default=0.95)
    ap.add_argument("--opacity_refine", type=float, default=0.95)

    ap.add_argument(
        "--show_noise",
        action="store_true",
        help="Show DBSCAN noise points (-1) in gray",
    )
    ap.add_argument("--noise_point_size", type=int, default=3)
    ap.add_argument("--noise_opacity", type=float, default=0.15)
    ap.add_argument("--noise_downsample", type=int, default=6)

    # Labels + summaries
    ap.add_argument(
        "--label_top", type=int, default=10, help="Annotate top-N clusters at centroid"
    )
    ap.add_argument("--label_font_size", type=int, default=14)
    ap.add_argument("--label_text_color", default="black")

    ap.add_argument(
        "--summary_top",
        type=int,
        default=10,
        help="Show top-N cluster ids/sizes in corner text",
    )
    ap.add_argument("--summary_font_size", type=int, default=12)

    # Save a screenshot (optional)
    ap.add_argument(
        "--screenshot",
        default=None,
        help="If set, save screenshot to this path instead of interactive show()",
    )

    args = ap.parse_args()

    runs_root = Path(args.runs_root)
    run_dir = (
        Path(args.run_dir)
        if args.run_dir
        else find_latest_run_dir(runs_root, args.rcsb)
    )
    print(f"[viz] Using run: {run_dir}")

    stage = run_dir / "stage"
    st10 = stage / "10_landmarks"
    st20 = stage / "20_exterior_shell"
    st40 = stage / "40_empty_space"
    st50 = stage / "50_clustering"
    st70 = stage / "70_mesh_validate"

    ptc, constr = load_landmarks(st10)

    alpha_shell = find_alpha_shell(st20)
    roi_json = find_roi_json(st40)
    final_mesh = find_final_mesh(st70)

    coarse_dir = st50 / "coarse"
    refine_dir = st50 / "refine"

    # Create a single window with 3 subplots
    pl = pv.Plotter(shape=(1, 3), border=True, window_size=(2100, 900))
    pl.set_background("white")

    # ---- Subplot 1: coarse
    pl.subplot(0, 0)
    pl.add_axes()
    pl.show_grid()
    if args.show_shell and alpha_shell is not None:
        add_mesh(pl, alpha_shell, opacity=args.shell_opacity, label="alpha_shell")
    elif args.show_shell:
        print("[viz] alpha shell not found; skipping")

    if args.show_roi:
        if roi_json is not None:
            add_roi_bbox(
                pl,
                roi_json,
                fallback_ptc=ptc,
                fallback_constr=constr,
                color=args.roi_color,
                line_width=args.roi_line_width,
                label="ROI",
            )
        else:
            print("[viz] ROI json not found; skipping ROI overlay")

    add_dbscan_pass_subplot(
        pl,
        coarse_dir,
        title="DBSCAN COARSE",
        point_size=args.point_size_coarse,
        opacity=args.opacity_coarse,
        downsample=max(1, args.downsample_coarse),
        show_noise=args.show_noise,
        noise_point_size=args.noise_point_size,
        noise_opacity=args.noise_opacity,
        noise_downsample=max(1, args.noise_downsample),
        render_spheres=args.render_spheres,
        label_top=args.label_top,
        label_font_size=args.label_font_size,
        label_text_color=args.label_text_color,
        summary_top=args.summary_top,
        summary_font_size=args.summary_font_size,
    )

    # ---- Subplot 2: refine
    pl.subplot(0, 1)
    pl.add_axes()
    pl.show_grid()
    if args.show_shell and alpha_shell is not None:
        add_mesh(pl, alpha_shell, opacity=args.shell_opacity, label="alpha_shell")

    if args.show_roi and roi_json is not None:
        add_roi_bbox(
            pl,
            roi_json,
            fallback_ptc=ptc,
            fallback_constr=constr,
            color=args.roi_color,
            line_width=args.roi_line_width,
            label="ROI",
        )

    add_dbscan_pass_subplot(
        pl,
        refine_dir,
        title="DBSCAN REFINE",
        point_size=args.point_size_refine,
        opacity=args.opacity_refine,
        downsample=max(1, args.downsample_refine),
        show_noise=args.show_noise,
        noise_point_size=args.noise_point_size,
        noise_opacity=args.noise_opacity,
        noise_downsample=max(1, args.noise_downsample),
        render_spheres=args.render_spheres,
        label_top=args.label_top,
        label_font_size=args.label_font_size,
        label_text_color=args.label_text_color,
        summary_top=args.summary_top,
        summary_font_size=args.summary_font_size,
    )

    # ---- Subplot 3: final (refine + mesh overlay)
    pl.subplot(0, 2)
    pl.add_axes()
    pl.show_grid()
    if args.show_shell and alpha_shell is not None:
        add_mesh(pl, alpha_shell, opacity=args.shell_opacity, label="alpha_shell")

    if args.show_roi and roi_json is not None:
        add_roi_bbox(
            pl,
            roi_json,
            fallback_ptc=ptc,
            fallback_constr=constr,
            color=args.roi_color,
            line_width=args.roi_line_width,
            label="ROI",
        )

    add_dbscan_pass_subplot(
        pl,
        refine_dir,
        title="FINAL (REFINE + MESH)",
        point_size=args.point_size_refine,
        opacity=min(1.0, args.opacity_refine),
        downsample=max(1, args.downsample_refine),
        show_noise=False,  # keep final clean by default
        noise_point_size=args.noise_point_size,
        noise_opacity=args.noise_opacity,
        noise_downsample=max(1, args.noise_downsample),
        render_spheres=args.render_spheres,
        label_top=args.label_top,
        label_font_size=args.label_font_size,
        label_text_color=args.label_text_color,
        summary_top=args.summary_top,
        summary_font_size=args.summary_font_size,
    )

    if args.show_mesh:
        if final_mesh is not None:
            add_mesh(pl, final_mesh, opacity=args.mesh_opacity, label="tunnel_mesh")
        else:
            print("[viz] final mesh not found; skipping")

    if args.link_views:
        try:
            pl.link_views()
        except Exception as e:
            print(f"[viz] link_views failed: {e}")

    # Camera / render
    pl.camera_position = "iso"

    if args.screenshot:
        out = Path(args.screenshot)
        out.parent.mkdir(parents=True, exist_ok=True)
        pl.show(screenshot=str(out), auto_close=True)
        print(f"[viz] wrote screenshot: {out}")
    else:
        pl.add_legend(bcolor="white", size=(0.18, 0.18))
        pl.show()


if __name__ == "__main__":
    main()
