import argparse
import json
import os
from pathlib import Path
import numpy as np
import pyvista as pv
import open3d as o3d
from ribctl.ribosome_ops import RibosomeOps

def visualize_artifact(artifact_path):
    """Visualize an artifact using the appropriate tool based on its file type"""
    path = Path(artifact_path)
    
    if not path.exists():
        print(f"File not found: {path}")
        return
    
    # Determine file type and visualize accordingly
    if path.suffix.lower() == '.npy':
        # Assume it's a point cloud saved as numpy array
        try:
            points = np.load(path)
            
            if len(points) == 0:
                print(f"Empty point cloud in {path}")
                return
                
            if points.ndim != 2 or points.shape[1] != 3:
                print(f"Warning: Expected 3D points, but found shape {points.shape}")
                if points.ndim == 1 and len(points) == 3:
                    # Single point case
                    points = points.reshape(1, 3)
                else:
                    print("Cannot visualize as point cloud")
                    return
            
            cloud = pv.PolyData(points)
            
            p = pv.Plotter()
            p.add_points(cloud, render_points_as_spheres=True, point_size=5)
            p.add_axes()
            p.show_grid()
            p.add_text(f"Point Cloud: {path.name} ({len(points)} points)", position="upper_edge", font_size=14)
            p.show()
        except Exception as e:
            print(f"Error visualizing NPY file: {e}")
        
    elif path.suffix.lower() == '.ply':
        try:
            if 'normal' in path.name.lower():
                # Assume it's a normal-estimated point cloud
                pcd = o3d.io.read_point_cloud(str(path))
                if len(pcd.points) == 0:
                    print(f"Empty point cloud in {path}")
                    return
                    
                # Check if normals exist
                if not pcd.has_normals():
                    print(f"Warning: Normal point cloud has no normals: {path}")
                    
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
        except Exception as e:
            print(f"Error visualizing PLY file: {e}")
    
    elif path.suffix.lower() in ['.png', '.jpg', '.jpeg']:
        # For images, just use the system's default viewer
        try:
            import subprocess
            import platform
            
            system = platform.system()
            if system == 'Darwin':  # macOS
                subprocess.run(['open', path])
            elif system == 'Windows':
                os.startfile(path)
            else:  # Linux
                subprocess.run(['xdg-open', path])
        except Exception as e:
            print(f"Error opening image file: {e}")
    
    else:
        print(f"Unsupported file type: {path.suffix}")

def find_structure_log(rcsb_id, log_dir=None):
    """
    Find the processing log for a structure, trying multiple potential locations.
    
    Args:
        rcsb_id: The structure ID
        log_dir: Optional override for log directory
        
    Returns:
        Path to the log file if found, None otherwise
    """
    # Normalize rcsb_id
    rcsb_id = rcsb_id.upper()
    
    potential_locations = []
    
    # First try the provided log_dir
    if log_dir:
        potential_locations.append(Path(log_dir) / f"{rcsb_id}_processing_log.json")
    
    # Try the structure-specific artifacts directory
    try:
        ro = RibosomeOps(rcsb_id)
        structure_dir = Path(ro.assets.paths.dir)
        
        # Structure's main artifacts directory
        potential_locations.append(structure_dir / "artifacts" / f"{rcsb_id}_processing_log.json")
        # Structure's directory
        potential_locations.append(structure_dir / f"{rcsb_id}_processing_log.json")
        # Default logs directory within structure dir
        potential_locations.append(structure_dir / "logs" / f"{rcsb_id}_processing_log.json")
    except Exception as e:
        print(f"Warning: Could not access structure directory: {e}")
    
    # Try default logs directory
    potential_locations.append(Path("/Users/rtviii/dev/riboxyz/ribctl/lib/npet/pipeline/logs") / f"{rcsb_id}_processing_log.json")
    
    # Check each location
    for loc in potential_locations:
        if loc.exists():
            print(f"Found log at: {loc}")
            return loc
    
    print(f"No processing log found for {rcsb_id}")
    print("Checked locations:")
    for loc in potential_locations:
        print(f"  - {loc}")
    
    return None

def visualize_structure_artifacts(rcsb_id, log_dir=None):
    """Visualize all artifacts for a structure"""
    log_path = find_structure_log(rcsb_id, log_dir)
    
    if not log_path:
        return
    
    # Load the log file
    try:
        with open(log_path, 'r') as f:
            log_data = json.load(f)
    except Exception as e:
        print(f"Error loading log file: {e}")
        return
    
    # Collect all artifacts
    artifacts = []
    artifact_count = 0
    
    for stage_name, stage_info in log_data["stages"].items():
        for artifact_path in stage_info.get("artifacts", []):
            artifact_count += 1
            if artifact_path and os.path.exists(artifact_path):
                artifacts.append({
                    "path": artifact_path,
                    "stage": stage_name,
                    "name": os.path.basename(artifact_path)
                })
            elif artifact_path:
                print(f"Warning: Artifact does not exist: {artifact_path}")
    
    if not artifacts:
        print(f"No accessible artifacts found for {rcsb_id}")
        print(f"Log references {artifact_count} artifacts that could not be found")
        return
    
    # Group artifacts by stage for better organization
    stages = {}
    for artifact in artifacts:
        stage = artifact["stage"]
        if stage not in stages:
            stages[stage] = []
        stages[stage].append(artifact)
    
    # Show a menu of artifacts organized by stage
    print(f"\nArtifacts for {rcsb_id}:")
    i = 1
    artifact_indices = {}
    
    for stage_name in sorted(stages.keys()):
        print(f"\n[{stage_name}]:")
        for artifact in stages[stage_name]:
            print(f"  {i}. {artifact['name']}")
            artifact_indices[i] = artifact
            i += 1
    
    # Let user select an artifact to visualize
    print("\nOptions:")
    print("  Enter a number to visualize that artifact")
    print("  Enter 's:[stage_name]' to visualize all artifacts in a stage")
    print("  Enter 'a' to visualize all artifacts sequentially")
    print("  Enter 'q' to quit")
    
    choice = input("\nYour choice: ")
    
    if choice.lower() == 'q':
        return
    elif choice.lower() == 'a':
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
    elif choice.lower().startswith('s:'):
        # Visualize all artifacts in a stage
        stage_name = choice[2:].strip()
        if stage_name in stages:
            print(f"\nVisualizing all artifacts in stage '{stage_name}':")
            for artifact in stages[stage_name]:
                print(f"\nVisualizing: {artifact['name']}")
                try:
                    visualize_artifact(artifact['path'])
                except Exception as e:
                    print(f"Error visualizing {artifact['path']}: {str(e)}")
                
                # Ask to continue after each visualization
                if stages[stage_name].index(artifact) < len(stages[stage_name]) - 1:
                    cont = input("Press Enter to continue to next artifact (or 'q' to quit): ")
                    if cont.lower() == 'q':
                        break
        else:
            print(f"Stage '{stage_name}' not found. Available stages: {', '.join(stages.keys())}")
    else:
        try:
            idx = int(choice)
            if idx in artifact_indices:
                artifact = artifact_indices[idx]
                print(f"\nVisualizing: {artifact['name']}")
                visualize_artifact(artifact['path'])
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
    vis_parser.add_argument("--log-dir", help="Directory containing logs (optional)")
    
    # Command to visualize a specific artifact
    file_parser = subparsers.add_parser("file", help="Visualize a specific artifact file")
    file_parser.add_argument("path", help="Path to the artifact file")
    
    # Command to list available structures
    list_parser = subparsers.add_parser("list", help="List structures with artifacts")
    list_parser.add_argument("--log-dir", help="Directory containing logs (optional)")
    
    args = parser.parse_args()
    
    if args.command == "structure":
        visualize_structure_artifacts(args.rcsb_id, args.log_dir)
    elif args.command == "file":
        visualize_artifact(args.path)
    elif args.command == "list":
        # You could add a function to list available structures with logs
        print("Not implemented yet")
    else:
        parser.print_help()

if __name__ == "__main__":
    main()