#!/usr/bin/env python3
"""
Dashboard and Artifact Viewer Integration

This script combines the dashboard and artifact viewer to create a seamless workflow:
1. Select structures in the dashboard
2. Visualize their artifacts in the gallery

Usage:
    python view_pipeline_artifacts.py
"""

import os
import sys
import webbrowser
import subprocess
from pathlib import Path

# Directory where this script is located
SCRIPT_DIR = Path(__file__).parent.absolute()

def run_dashboard():
    """Run the dashboard generator and open it in a browser"""
    print("Generating dashboard...")
    
    try:
        # First try to find dashboard script in local directory
        dashboard_script = SCRIPT_DIR / "html_tally.py"
        
        if not dashboard_script.exists():
            # Try in the lib/npet/pipeline directory
            pipeline_dir = Path("/Users/rtviii/dev/riboxyz/ribctl/lib/npet/pipeline")
            dashboard_script = pipeline_dir / "html_tally.py"
        
        if not dashboard_script.exists():
            print("Could not find dashboard script (html_tally.py)")
            print("Please provide the path to the dashboard script:")
            path_input = input()
            dashboard_script = Path(path_input)
            
            if not dashboard_script.exists():
                raise FileNotFoundError(f"Dashboard script not found at {dashboard_script}")
        
        # Run dashboard generator
        log_dir = Path("/Users/rtviii/dev/riboxyz/ribctl/lib/npet/pipeline/logs")
        output_file = SCRIPT_DIR / "dashboard.html"
        
        # Ask user if they want to specify a different log directory
        print(f"Using log directory: {log_dir}")
        print("Enter different log directory path or press Enter to continue:")
        log_dir_input = input()
        if log_dir_input:
            log_dir = Path(log_dir_input)
        
        # Ensure log directory exists
        if not log_dir.exists():
            raise FileNotFoundError(f"Log directory not found: {log_dir}")
        
        # Run dashboard generator
        cmd = [sys.executable, str(dashboard_script), str(log_dir), str(output_file)]
        print(f"Running: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print("Error running dashboard generator:")
            print(result.stderr)
            return None
        
        # Open dashboard in browser
        dashboard_url = f"file://{output_file.absolute()}"
        print(f"Opening dashboard: {dashboard_url}")
        webbrowser.open(dashboard_url)
        
        return output_file
    
    except Exception as e:
        print(f"Error generating dashboard: {e}")
        return None

def run_artifact_viewer(structures, output_dir, artifact_type_args=None):
    """Run the artifact viewer for selected structures"""
    print("Generating artifact gallery...")
    
    try:
        # Find artifact viewer script
        viewer_script = SCRIPT_DIR / "bulk_artifact_viewer.py"
        
        if not viewer_script.exists():
            # Try in the current directory
            viewer_script = Path("bulk_artifact_viewer.py")
        
        if not viewer_script.exists():
            print("Could not find artifact viewer script (bulk_artifact_viewer.py)")
            print("Please provide the path to the artifact viewer script:")
            path_input = input()
            viewer_script = Path(path_input)
            
            if not viewer_script.exists():
                raise FileNotFoundError(f"Artifact viewer script not found at {viewer_script}")
        
        # Prepare command
        cmd = [
            sys.executable, 
            str(viewer_script), 
            "--structures"
        ] + structures + [
            "--gallery",
            "--parallel",
            "--output-dir", output_dir
        ]
        
        # Add artifact type args if specified
        if artifact_type_args:
            cmd.extend(artifact_type_args)
            
        print(f"Running: {' '.join(cmd)}")
        
        # Run viewer
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        if result.returncode != 0:
            print("Error running artifact viewer")
            return False
        
        return True
    
    except Exception as e:
        print(f"Error running artifact viewer: {e}")
        return False

def main():
    """Main function"""
    print("=== Pipeline Artifact Explorer ===")
    print("This tool helps you select structures and visualize their artifacts.\n")
    
    # Run dashboard
    dashboard_file = run_dashboard()
    
    if not dashboard_file:
        print("Could not generate dashboard. Exiting.")
        return
    
    # Wait for user to select structures
    print("\nInstructions:")
    print("1. In the dashboard, click on rows to select structures")
    print("2. Copy the list of selected structures from the panel at the bottom")
    print("3. Return to this window and paste the list when prompted")
    print("\nPress Enter when ready to continue...")
    input()
    
    # Get selected structures
    print("\nPaste the structure IDs you selected:")
    structures_input = input()
    
    # Parse structures (handle various formats)
    structures = []
    for line in structures_input.splitlines():
        structures.extend(line.strip().split())
    
    if not structures:
        print("No structures provided. Exiting.")
        return
    
    print(f"Selected {len(structures)} structures: {', '.join(structures[:5])}" + 
          (f" and {len(structures)-5} more..." if len(structures) > 5 else ""))
    
    # Ask which artifact types to visualize
    print("\nWhich artifacts would you like to visualize?")
    print("1. All artifacts")
    print("2. Only mesh files (.ply)")
    print("3. Only visualizations (.png)")
    print("4. Only point clouds (.npy)")
    
    artifact_choice = input("Enter your choice (1-4): ")
    
    artifact_type_arg = []
    if artifact_choice == "2":
        artifact_type_arg = ["--type", ".ply"]
    elif artifact_choice == "3":
        artifact_type_arg = ["--type", ".png"]
    elif artifact_choice == "4":
        artifact_type_arg = ["--type", ".npy"]
    
    # Ask for output directory
    print("\nWhere would you like to save the gallery? (default: ./pipeline_gallery)")
    output_dir = input("Enter path or press Enter for default: ")
    
    if not output_dir:
        output_dir = "./pipeline_gallery"
    
    # Ensure output directory exists
    Path(output_dir).mkdir(exist_ok=True, parents=True)
    
    # Create structure and stage subdirectories
    viewer_script = SCRIPT_DIR / "bulk_artifact_viewer.py"
    
    if not viewer_script.exists():
        # Try in the current directory
        viewer_script = Path("bulk_artifact_viewer.py")
    
    if not viewer_script.exists():
        print("Could not find artifact viewer script (bulk_artifact_viewer.py)")
        print("Please provide the path to the artifact viewer script:")
        path_input = input()
        viewer_script = Path(path_input)
        
        if not viewer_script.exists():
            raise FileNotFoundError(f"Artifact viewer script not found at {viewer_script}")
    
    # Run artifact viewer with specified output directory
    cmd = [
        sys.executable, 
        str(viewer_script), 
        "--structures"
    ] + structures + [
        "--gallery",
        "--parallel",
        "--output-dir", output_dir
    ] + artifact_type_arg
    
    print(f"Running: {' '.join(cmd)}")
    
    try:
        # Run viewer
        result = subprocess.run(cmd, capture_output=False, text=True)
        
        if result.returncode != 0:
            print("\nError generating artifact gallery.")
        else:
            print("\nArtifact gallery generated successfully!")
            print(f"Gallery saved to: {os.path.abspath(output_dir)}/gallery.html")
            
    except Exception as e:
        print(f"\nError running artifact viewer: {e}")

if __name__ == "__main__":
    main()