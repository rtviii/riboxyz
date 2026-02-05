# __scripts/npet2_slice_viewer.py
#!/usr/bin/env python3
"""
Slice viewer for NPET2 voxel grids.
Press 'x', 'y', or 'z' to change axis. Arrow keys to move slice.
"""

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def find_latest_run_dir(runs_root: Path, rcsb: str) -> Path:
    rdir = runs_root / rcsb.upper()
    candidates = [p for p in rdir.iterdir() if p.is_dir()]
    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0]


def load_grid(base_path: Path):
    data_path = base_path.parent / f"{base_path.stem}_data.npy"
    spec_path = base_path.parent / f"{base_path.stem}_spec.json"
    
    data = np.load(data_path)
    spec = json.loads(spec_path.read_text())
    return spec, data


class SliceViewer:
    def __init__(self, data: np.ndarray, title: str = "Grid Slice"):
        self.data = data
        self.axis = 0  # 0=X, 1=Y, 2=Z
        self.index = data.shape[self.axis] // 2
        
        self.fig, self.ax = plt.subplots(figsize=(10, 10))
        self.fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.title = title
        
        self.update()
    
    def get_slice(self):
        if self.axis == 0:
            return self.data[self.index, :, :]
        elif self.axis == 1:
            return self.data[:, self.index, :]
        else:
            return self.data[:, :, self.index]
    
    def update(self):
        self.ax.clear()
        slice_data = self.get_slice()
        
        self.ax.imshow(slice_data.T, origin='lower', cmap='gray', interpolation='nearest')
        
        axis_names = ['X', 'Y', 'Z']
        self.ax.set_title(f"{self.title} | {axis_names[self.axis]} = {self.index}/{self.data.shape[self.axis]-1}")
        self.ax.set_xlabel('Voxel index')
        self.ax.set_ylabel('Voxel index')
        
        # Add occupancy percentage
        occ_pct = 100.0 * slice_data.sum() / slice_data.size
        self.ax.text(0.02, 0.98, f"Slice occupancy: {occ_pct:.1f}%", 
                    transform=self.ax.transAxes, va='top', 
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        
        self.fig.canvas.draw()
    
    def on_key(self, event):
        if event.key == 'x':
            self.axis = 0
            self.index = self.data.shape[self.axis] // 2
        elif event.key == 'y':
            self.axis = 1
            self.index = self.data.shape[self.axis] // 2
        elif event.key == 'z':
            self.axis = 2
            self.index = self.data.shape[self.axis] // 2
        elif event.key == 'left':
            self.index = max(0, self.index - 1)
        elif event.key == 'right':
            self.index = min(self.data.shape[self.axis] - 1, self.index + 1)
        elif event.key == 'down':
            self.index = max(0, self.index - 10)
        elif event.key == 'up':
            self.index = min(self.data.shape[self.axis] - 1, self.index + 10)
        else:
            return
        
        self.update()


def main():
    ap = argparse.ArgumentParser(description="Slice viewer for NPET2 grids")
    ap.add_argument("--runs_root", default="/Users/rtviii/dev/riboxyz/NPET2/runs")
    ap.add_argument("--rcsb", required=True)
    ap.add_argument("--run_dir", default=None)
    ap.add_argument("--grid", default="occupancy_grid_level_0", 
                   help="Grid name (e.g., occupancy_grid_level_0, empty_mask_level_1)")
    
    args = ap.parse_args()
    
    runs_root = Path(args.runs_root)
    if args.run_dir:
        run_dir = Path(args.run_dir)
    else:
        run_dir = find_latest_run_dir(runs_root, args.rcsb)
    
    grid_path = run_dir / "stage" / "40_empty_space" / args.grid
    
    print(f"Loading grid: {grid_path}")
    spec, data = load_grid(grid_path)
    
    print(f"Grid shape: {data.shape}")
    print(f"Voxel size: {spec['voxel_size']} Å")
    print(f"Occupancy: {100.0 * data.sum() / data.size:.1f}%")
    print("\nControls:")
    print("  x/y/z: Switch axis")
    print("  Left/Right: Move slice ±1")
    print("  Up/Down: Move slice ±10")
    
    viewer = SliceViewer(data, title=f"{args.rcsb} - {args.grid}")
    plt.show()


if __name__ == "__main__":
    main()