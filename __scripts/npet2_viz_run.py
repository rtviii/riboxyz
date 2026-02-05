import argparse
import numpy as np
import pyvista as pv
from pathlib import Path

def load_npy(p): return np.load(p)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    args = ap.parse_args()
    run = Path(args.run_dir)

    shell = pv.read(run / "stage/20_exterior_shell/alpha_shell.ply")
    tunnel = pv.read(run / "stage/70_mesh_validate/npet2_tunnel_mesh.ply")

    p40 = run / "stage/40_empty_space"
    npys = sorted(p40.glob("empty_points_*.npy"))

    pl = pv.Plotter()
    pl.add_mesh(shell, opacity=0.15)
    pl.add_mesh(tunnel, opacity=0.45)

    for p in npys:
        pts = load_npy(p)
        cloud = pv.PolyData(pts)
        pl.add_mesh(cloud, point_size=3, render_points_as_spheres=True, label=p.name)

    pl.add_legend()
    pl.show()

if __name__ == "__main__":
    main()
