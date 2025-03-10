import argparse
import json
import os
from pathlib import Path
import numpy as np
import pyvista as pv
import open3d as o3d

def visualize_artifact(artifact_path):
    """Visualize an artifact using the appropriate tool based on its file type"""
    path = Path(artifact_path)
    
    if not path.exists():
        print(f"File not found: {path}")
        return
    
    # Determine file type and visualize accordingly
    if path.suffix.lower() == '.npy':
        # Assume it's a point cloud saved as numpy array
        points = np.load(path)
        cloud = pv.PolyData(points)
        
        p = pv.Plotter()
        p.add_points(cloud, render_points_as_spheres=True, point_size=5)
        p.add_axes()
        p.show_grid()
        p.add_text(f"Point Cloud: {path.name}", position="upper_edge", font_size=14)
        p.show()
        
    elif path.suffix.lower() == '.ply':
        if 'normal' in path.name.lower():
            # Assume it's a normal-estimated point cloud
            pcd = o3d.io.read_point_cloud(str(path))
            o3d.visualization.draw_geometries([pcd], point_show_normal=True)
        else:
            # Assume it's a mesh
            mesh = pv.read(path)
            p = pv.Plotter()
            p.add_mesh(mesh, show_edges=True)
            p.add_axes()
            p.show_grid()
            p.add_text(f"Mesh: {path.name}", position="upper_edge", font_size=14)
            p.show()
    
    elif path.suffix.lower() in ['.png', '.jpg', '.jpeg']:
        # For images, just use the system's default viewer
        import subprocess
        import platform
        
        system = platform.system()
        if system == 'Darwin':  # macOS
            subprocess.run(['open', path])
        elif system == 'Windows':
            os.startfile(path)
        else:  # Linux
            subprocess.run(['xdg-open', path])
    
    else:
        print(f"Unsupported file type: {path.suffix}")

def visualize_structure_artifacts(rcsb_id, log_dir):
    """Visualize all artifacts for a structure"""
    log_path = Path(log_dir) / f"{rcsb_id.upper()}_processing_log.json"
    
    if not log_path.exists():
        print(f"No processing log found for {rcsb_id} at {log_path}")
        return
    
    # Load the log file
    with open(log_path, 'r') as f:
        log_data = json.load(f)
    
    # Collect all artifacts
    artifacts = []
    for stage_name, stage_info in log_data["stages"].items():
        for artifact_path in stage_info.get("artifacts", []):
            if artifact_path and os.path.exists(artifact_path):
                artifacts.append({
                    "path": artifact_path,
                    "stage": stage_name,
                    "name": os.path.basename(artifact_path)
                })
    
    if not artifacts:
        print(f"No artifacts found for {rcsb_id}")
        return
    
    # Show a menu of artifacts
    print(f"\nArtifacts for {rcsb_id}:")
    for i, artifact in enumerate(artifacts, 1):
        print(f"{i}. [{artifact['stage']}] {artifact['name']}")
    
    # Let user select an artifact to visualize
    choice = input("\nEnter number to visualize (or 'a' for all): ")
    
    if choice.lower() == 'a':
        # Visualize all artifacts one by one
        for artifact in artifacts:
            print(f"\nVisualizing: {artifact['name']}")
            try:
                visualize_artifact(artifact['path'])
            except Exception as e:
                print(f"Error visualizing {artifact['path']}: {str(e)}")
            
            # Ask to continue after each visualization
            if artifacts.index(artifact) < len(artifacts) - 1:
                cont = input("Press Enter to continue to next artifact (or 'q' to quit): ")
                if cont.lower() == 'q':
                    break
    else:
        try:
            idx = int(choice) - 1
            if 0 <= idx < len(artifacts):
                visualize_artifact(artifacts[idx]['path'])
            else:
                print("Invalid selection")
        except ValueError:
            print("Invalid input")

def main():
    parser = argparse.ArgumentParser(description="NPET Artifact Visualization Tool")
    
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # Command to visualize artifacts for a structure
    vis_parser = subparsers.add_parser("structure", help="Visualize artifacts for a structure")
    vis_parser.add_argument("rcsb_id", help="RCSB ID of the structure")
    vis_parser.add_argument("--log-dir", default="logs", help="Directory containing logs")
    
    # Command to visualize a specific artifact
    file_parser = subparsers.add_parser("file", help="Visualize a specific artifact file")
    file_parser.add_argument("path", help="Path to the artifact file")
    
    args = parser.parse_args()
    
    if args.command == "structure":
        visualize_structure_artifacts(args.rcsb_id, args.log_dir)
    elif args.command == "file":
        visualize_artifact(args.path)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()