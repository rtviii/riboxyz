#!/usr/bin/env python3
import os
import argparse
import numpy as np
import pyvista as pv
import open3d as o3d
from pathlib import Path
from ribctl.asset_manager.asset_types import AssetType
from ribctl.lib.npet.various_visualization import (
    visualize_filtered_residues,
    visualize_mesh,
    visualize_pointcloud,
    visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs,
)
from ribctl.lib.npet.tunnel_asset_manager import TunnelMeshAssetsManager
from ribctl.lib.npet.kdtree_approach import (
    landmark_ptc,
    landmark_constriction_site,
)

class NPETVisualizer:
    """
    Visualization tool for NPET tunnel artifacts without recomputing them.
    Allows inspection of successfully created tunnels and intermediate results.
    """
    
    def __init__(self, rcsb_id):
        """Initialize with RCSB ID for the structure"""
        self.rcsb_id = rcsb_id
        self.assets = TunnelMeshAssetsManager(rcsb_id)
        
        # Key landmark points
        self.ptc_pt = np.array(landmark_ptc(rcsb_id))
        self.constriction_pt = np.array(landmark_constriction_site(rcsb_id))
        
        # Asset paths
        self.cifpath = AssetType.MMCIF.get_path(rcsb_id)
        self.ashapepath = AssetType.ALPHA_SHAPE.get_path(rcsb_id)
        self.meshpath = AssetType.NPET_MESH.get_path(rcsb_id)
        
        # Pipeline parameters (from original function)
        self.R = 35
        self.H = 120
    
    def check_assets_exist(self):
        """Check if the required assets exist and return a status report"""
        assets_status = {
            "MMCIF": os.path.exists(self.cifpath),
            "Alpha Shape": os.path.exists(self.ashapepath),
            "NPET Mesh": os.path.exists(self.meshpath),
            "PCD Normal": os.path.exists(self.assets.tunnel_pcd_normal_estimated),
            "Filtered Residues": os.path.exists(self.assets.tunnel_filtered_residues),
            "DBSCAN Initial": os.path.exists(self.assets.tunnel_dbscan_initial),
            "DBSCAN Refined": os.path.exists(self.assets.tunnel_dbscan_refined),
            "Surface Points": os.path.exists(self.assets.tunnel_surface_points),
        }
        return assets_status
    
    def visualize_all(self):
        """Run all available visualizations"""
        status = self.check_assets_exist()
        print(f"Asset status for {self.rcsb_id}:")
        for asset, exists in status.items():
            print(f"  {asset}: {'✓' if exists else '✗'}")
        
        # Visualize each step that has data
        if status["NPET Mesh"]:
            self.visualize_final_mesh()
        
        if status["Surface Points"]:
            self.visualize_surface_points()
        
        if status["DBSCAN Refined"]:
            self.visualize_dbscan_refined()
        
        if status["DBSCAN Initial"]:
            self.visualize_dbscan_initial()
        
        if status["Filtered Residues"] and status["MMCIF"]:
            self.visualize_residues()
    
    def visualize_final_mesh(self):
        """Visualize the final NPET mesh"""
        print(f"\nVisualizing final mesh for {self.rcsb_id}")
        if os.path.exists(self.meshpath):
            visualize_mesh(self.meshpath, self.rcsb_id)
        else:
            print(f"Final mesh not found at {self.meshpath}")
    
    def visualize_surface_points(self):
        """Visualize surface points"""
        print(f"\nVisualizing surface points for {self.rcsb_id}")
        if os.path.exists(self.assets.tunnel_surface_points):
            points = np.load(self.assets.tunnel_surface_points)
            visualize_pointcloud(points, self.rcsb_id)
        else:
            print(f"Surface points not found at {self.assets.tunnel_surface_points}")
    
    def visualize_dbscan_refined(self):
        """Visualize refined DBSCAN clusters"""
        print(f"\nVisualizing refined DBSCAN clusters for {self.rcsb_id}")
        if os.path.exists(self.assets.tunnel_dbscan_refined):
            refined_cluster = np.load(self.assets.tunnel_dbscan_refined)
            # Use 3.5 and 175 from the original script
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                {0: refined_cluster},  # Mock clusters container
                3.5, 175,
                self.ptc_pt, self.constriction_pt,
                refined_cluster, self.R, self.H
            )
        else:
            print(f"Refined DBSCAN clusters not found at {self.assets.tunnel_dbscan_refined}")
    
    def visualize_dbscan_initial(self):
        """Visualize initial DBSCAN clusters"""
        print(f"\nVisualizing initial DBSCAN clusters for {self.rcsb_id}")
        if os.path.exists(self.assets.tunnel_dbscan_initial):
            largest_cluster = np.load(self.assets.tunnel_dbscan_initial)
            # Use 5.5 and 600 from the original script
            visualize_DBSCAN_CLUSTERS_particular_eps_minnbrs(
                {0: largest_cluster},  # Mock clusters container
                5.5, 600,
                self.ptc_pt, self.constriction_pt,
                largest_cluster, self.R, self.H
            )
        else:
            print(f"Initial DBSCAN clusters not found at {self.assets.tunnel_dbscan_initial}")
    
    def visualize_residues(self):
        """Visualize filtered residues"""
        print(f"\nVisualizing filtered residues for {self.rcsb_id}")
        if os.path.exists(self.assets.tunnel_filtered_residues):
            # Would need to load residues data properly
            # This is simplified and may need adaptation based on how you store residues
            print("Note: Residue visualization requires loading Bio.PDB structures")
            print("You may need to modify this function based on how residues are stored")
        else:
            print(f"Filtered residues not found at {self.assets.tunnel_filtered_residues}")
    
    def visualize_open3d_point_cloud(self):
        """Visualize Open3D point cloud with normals if available"""
        print(f"\nVisualizing point cloud with normals for {self.rcsb_id}")
        if os.path.exists(self.assets.tunnel_pcd_normal_estimated):
            pcd = o3d.io.read_point_cloud(self.assets.tunnel_pcd_normal_estimated)
            o3d.visualization.draw_geometries([pcd], 
                                             window_name=f"NPET Point Cloud - {self.rcsb_id}",
                                             point_show_normal=True)
        else:
            print(f"Point cloud with normals not found at {self.assets.tunnel_pcd_normal_estimated}")

def check_all_structures(base_dir=None):
    """Check all structures in a directory for which we have computed NPET meshes"""
    if base_dir is None:
        # Try to get base directory from environment or use default
        base_dir = os.environ.get("RIBCTL_DATA_DIR", ".")
    
    # Find all NPET mesh files
    mesh_dir = os.path.join(base_dir, "structures")
    successful_structures = []
    
    for root, dirs, files in os.walk(mesh_dir):
        for file in files:
            if file.endswith("_npet_mesh.ply"):
                rcsb_id = file.split("_npet_mesh.ply")[0]
                successful_structures.append(rcsb_id)
    
    print(f"Found {len(successful_structures)} structures with NPET meshes")
    return successful_structures


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="NPET Tunnel Artifact Visualization Tool")
    parser.add_argument("--rcsb_id", type=str, help="RCSB ID to visualize")
    parser.add_argument("--list", action="store_true", help="List all structures with NPET meshes")
    parser.add_argument("--all", action="store_true", help="Visualize all structures with NPET meshes")
    parser.add_argument("--step", choices=["mesh", "points", "dbscan_refined", "dbscan_initial", "residues", "pcd"],
                        help="Visualize specific step only")
    
    args = parser.parse_args()
    
    if args.list:
        structures = check_all_structures()
        print("Structures with NPET meshes:")
        for struct in structures:
            print(f"  {struct}")
        exit(0)
    
    if args.all:
        structures = check_all_structures()
        for struct in structures:
            print(f"\n{'='*50}")
            print(f"Visualizing structure: {struct}")
            print(f"{'='*50}")
            visualizer = NPETVisualizer(struct)
            visualizer.visualize_all()
        exit(0)
    
    if not args.rcsb_id:
        parser.print_help()
        exit(1)
    
    visualizer = NPETVisualizer(args.rcsb_id)
    
    if args.step:
        # Call specific visualization function
        if args.step == "mesh":
            visualizer.visualize_final_mesh()
        elif args.step == "points":
            visualizer.visualize_surface_points()
        elif args.step == "dbscan_refined":
            visualizer.visualize_dbscan_refined()
        elif args.step == "dbscan_initial":
            visualizer.visualize_dbscan_initial()
        elif args.step == "residues":
            visualizer.visualize_residues()
        elif args.step == "pcd":
            visualizer.visualize_open3d_point_cloud()
    else:
        # Run all visualizations
        visualizer.visualize_all()