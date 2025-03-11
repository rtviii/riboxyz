#!/usr/bin/env python3
"""
Bulk Artifact Visualization Tool

This script allows visualizing artifacts from multiple structures at once.
It integrates with the pipeline dashboard to select structures and view their artifacts.

Usage:
    python bulk_viewer.py [--log-dir LOG_DIR] [--output-dir OUTPUT_DIR] [--max-structures MAX]
"""

import os
import sys
import json
import argparse
import tempfile
from pathlib import Path
import subprocess
import concurrent.futures
import numpy as np
import pyvista as pv
import open3d as o3d
from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import shutil
import webbrowser
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor

try:
    from ribctl.ribosome_ops import RibosomeOps
except ImportError:
    print("Warning: ribctl not available, some features may be limited")
    RibosomeOps = None

# Default artifact types for each stage that are interesting to visualize
DEFAULT_VISUALIZATION_ARTIFACTS = {
    "alpha_shape": ["_ALPHA_SHAPE.ply", "_alpha_mesh.png"],
    "entity_filtering": ["_filtered_residues.png"],
    "point_cloud_processing": ["_point_cloud.png", "_interior_points.npy"],
    "clustering": ["_clusters.png", "_largest_cluster.npy"],
    "refinement": ["_refined_clusters.png", "_refined_cluster.npy"],
    "surface_extraction": ["_surface_points.png", "_surface_points.npy"],
    "normal_estimation": ["_normal_estimated_pcd.ply"],
    "mesh_reconstruction": ["_NPET_MESH.ply"],
    "validation": ["_mesh.png"]
}

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="Bulk Artifact Visualization Tool")
    parser.add_argument("--log-dir", help="Directory containing logs")
    parser.add_argument("--output-dir", default="pipeline_gallery", help="Directory to store generated visualizations")
    parser.add_argument("--dashboard", action="store_true", help="Launch dashboard first to select structures")
    parser.add_argument("--structures", nargs="+", help="List of structure IDs to visualize")
    parser.add_argument("--stage", help="Specific stage to visualize artifacts from")
    parser.add_argument("--max-structures", type=int, default=20, help="Maximum number of structures to visualize")
    parser.add_argument("--gallery", action="store_true", help="Create HTML gallery of artifacts")
    parser.add_argument("--parallel", action="store_true", help="Process artifacts in parallel")
    parser.add_argument("--type", help="Filter artifacts by file type (e.g., .ply, .npy, .png)")
    
    return parser.parse_args()

def find_logs(log_dir):
    """Find all processing log files in the log directory"""
    log_dir = Path(log_dir)
    log_files = list(log_dir.glob("*_processing_log.json"))
    
    if not log_files:
        print(f"No log files found in {log_dir}")
        return []
    
    return log_files

def load_log(log_path):
    """Load and parse a log file"""
    try:
        with open(log_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading {log_path}: {e}")
        return None

def find_structure_logs(rcsb_ids, log_dir=None):
    """Find logs for specific structures"""
    logs = {}
    
    # Normalize IDs
    rcsb_ids = [id.upper() for id in rcsb_ids]
    
    # Potential log directories
    potential_dirs = []
    if log_dir:
        potential_dirs.append(Path(log_dir))
    
    # Default pipeline logs directory
    potential_dirs.append(Path("/Users/rtviii/dev/riboxyz/ribctl/lib/npet/pipeline/logs"))
    
    # Current directory and its logs subdirectory
    potential_dirs.append(Path.cwd())
    potential_dirs.append(Path.cwd() / "logs")
    
    # Try to find logs for each structure
    for rcsb_id in rcsb_ids:
        found = False
        
        # Try each potential directory
        for dir_path in potential_dirs:
            if not dir_path.exists():
                continue
                
            log_path = dir_path / f"{rcsb_id}_processing_log.json"
            if log_path.exists():
                log_data = load_log(log_path)
                if log_data:
                    logs[rcsb_id] = {
                        "log_path": log_path,
                        "data": log_data
                    }
                    found = True
                    break
        
        if not found and RibosomeOps:
            # Try to find logs in the structure's own directory
            try:
                ro = RibosomeOps(rcsb_id)
                structure_dir = Path(ro.assets.paths.dir)
                
                # Check in structure directory
                log_path = structure_dir / f"{rcsb_id}_processing_log.json"
                if log_path.exists():
                    log_data = load_log(log_path)
                    if log_data:
                        logs[rcsb_id] = {
                            "log_path": log_path,
                            "data": log_data
                        }
                        found = True
                
                if not found:
                    # Check in artifacts subdirectory
                    log_path = structure_dir / "artifacts" / f"{rcsb_id}_processing_log.json"
                    if log_path.exists():
                        log_data = load_log(log_path)
                        if log_data:
                            logs[rcsb_id] = {
                                "log_path": log_path,
                                "data": log_data
                            }
                            found = True
            except Exception as e:
                print(f"Could not access structure directory for {rcsb_id}: {e}")
        
        if not found:
            print(f"No log found for {rcsb_id}")
    
    return logs

def collect_artifacts(logs, stage=None):
    """
    Collect artifacts from logs, optionally filtering by stage
    
    Returns dict:
    {
        rcsb_id: {
            stage_name: [
                {path: artifact_path, exists: bool}
            ]
        }
    }
    """
    artifacts = {}
    
    for rcsb_id, log_info in logs.items():
        log_data = log_info["data"]
        artifacts[rcsb_id] = {}
        
        for stage_name, stage_data in log_data["stages"].items():
            if stage and stage != stage_name:
                continue
                
            if "artifacts" in stage_data and stage_data["artifacts"]:
                artifacts[rcsb_id][stage_name] = []
                
                for artifact_path in stage_data["artifacts"]:
                    # Skip empty paths
                    if not artifact_path:
                        continue
                        
                    path = Path(artifact_path)
                    exists = path.exists()
                    
                    # Store information about the artifact
                    artifacts[rcsb_id][stage_name].append({
                        "path": str(path),
                        "name": path.name,
                        "exists": exists,
                        "type": path.suffix.lower(),
                        "stage": stage_name,
                        "rcsb_id": rcsb_id
                    })
    
    return artifacts

def filter_artifacts_by_type(artifacts, artifact_types=None, file_type=None):
    """Filter artifacts to include only specific types/names"""
    if not artifact_types:
        # Default to showing the most interesting visualizable artifacts
        artifact_types = []
        for stage, patterns in DEFAULT_VISUALIZATION_ARTIFACTS.items():
            artifact_types.extend(patterns)
    
    filtered = {}
    for rcsb_id, stages in artifacts.items():
        filtered[rcsb_id] = {}
        
        for stage_name, stage_artifacts in stages.items():
            filtered_artifacts = []
            
            for artifact in stage_artifacts:
                # Apply pattern filter
                pattern_match = any(pattern in artifact["name"] for pattern in artifact_types)
                
                # Apply file type filter if specified
                type_match = True
                if file_type:
                    type_match = artifact["type"].endswith(file_type)
                
                if pattern_match and type_match:
                    filtered_artifacts.append(artifact)
                    
            if filtered_artifacts:
                filtered[rcsb_id][stage_name] = filtered_artifacts
    
    return filtered

def render_mesh_to_image(mesh_path, output_path=None, resolution=(800, 600), return_image=False):
    """Render a 3D mesh to a PNG image"""
    try:
        # Create PyVista plotter
        p = pv.Plotter(off_screen=True)
        p.window_size = resolution
        
        # Load and add mesh
        mesh = pv.read(mesh_path)
        p.add_mesh(mesh, show_edges=False, color="lightblue")
        p.add_axes()
        p.view_isometric()
        
        # Render and save
        if output_path:
            p.screenshot(output_path)
        
        if return_image:
            image_bytes = p.screenshot(return_img=True)
            return Image.fromarray(image_bytes)
        
        return True
        
    except Exception as e:
        print(f"Error rendering mesh {mesh_path}: {e}")
        return False

def render_point_cloud_to_image(points_path, output_path=None, resolution=(800, 600), return_image=False):
    """Render a point cloud to a PNG image"""
    try:
        # Load points
        points = np.load(points_path)
        
        # Check point cloud validity
        if len(points) == 0:
            print(f"Empty point cloud: {points_path}")
            return False
            
        if points.ndim != 2 or points.shape[1] != 3:
            print(f"Invalid point cloud shape {points.shape}: {points_path}")
            if points.ndim == 1 and len(points) == 3:
                # Single point case
                points = points.reshape(1, 3)
            else:
                return False
        
        # Create PyVista plotter
        p = pv.Plotter(off_screen=True)
        p.window_size = resolution
        
        # Add points
        cloud = pv.PolyData(points)
        p.add_points(cloud, render_points_as_spheres=True, point_size=5, color="lightblue")
        p.add_axes()
        p.view_isometric()
        
        # Render and save
        if output_path:
            p.screenshot(output_path)
        
        if return_image:
            image_bytes = p.screenshot(return_img=True)
            return Image.fromarray(image_bytes)
        
        return True
        
    except Exception as e:
        print(f"Error rendering point cloud {points_path}: {e}")
        return False

def render_normal_cloud_to_image(pcd_path, output_path=None, resolution=(800, 600), return_image=False):
    """Render a normal-estimated point cloud to a PNG image"""
    try:
        # Load with Open3D
        pcd = o3d.io.read_point_cloud(str(pcd_path))
        
        if len(pcd.points) == 0:
            print(f"Empty point cloud: {pcd_path}")
            return False
        
        # Convert to numpy arrays
        points = np.asarray(pcd.points)
        
        # Create PyVista plotter (O3D visualization can't render to image easily)
        p = pv.Plotter(off_screen=True)
        p.window_size = resolution
        
        # Add points
        cloud = pv.PolyData(points)
        p.add_points(cloud, render_points_as_spheres=True, point_size=5, color="lightblue")
        
        # Add a subset of normals as arrows if they exist
        if pcd.has_normals():
            normals = np.asarray(pcd.normals)
            # Sample a subset of points to show normals (avoid clutter)
            n_samples = min(500, len(points))
            indices = np.random.choice(len(points), n_samples, replace=False)
            
            sampled_points = points[indices]
            sampled_normals = normals[indices]
            
            # Scale normals for visibility
            normal_length = np.mean(np.linalg.norm(points, axis=1)) * 0.05
            scaled_normals = sampled_normals * normal_length
            
            p.add_arrows(sampled_points, scaled_normals, color="red")
        
        p.add_axes()
        p.view_isometric()
        
        # Render and save
        if output_path:
            p.screenshot(output_path)
        
        if return_image:
            image_bytes = p.screenshot(return_img=True)
            return Image.fromarray(image_bytes)
        
        return True
        
    except Exception as e:
        print(f"Error rendering normal point cloud {pcd_path}: {e}")
        return False

def process_artifact(artifact, output_dir=None):
    """
    Process a single artifact - either copy existing image or render 3D data to image
    
    Returns a path to the image representation of the artifact
    """
    artifact_path = artifact["path"]
    artifact_type = artifact["type"]
    
    # Skip if artifact doesn't exist
    if not os.path.exists(artifact_path):
        print(f"Artifact doesn't exist: {artifact_path}")
        return None
    
    # Determine output path if needed
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
        
        rcsb_id = artifact["rcsb_id"]
        stage = artifact["stage"]
        
        # Create structure and stage subdirectories
        struct_dir = output_dir / rcsb_id
        struct_dir.mkdir(exist_ok=True)
        
        stage_dir = struct_dir / stage
        stage_dir.mkdir(exist_ok=True)
        
        basename = Path(artifact_path).stem
        output_path = stage_dir / f"{basename}.png"
    else:
        # Use temporary file
        output_path = Path(tempfile.mktemp(suffix=".png"))
    
    # Process based on artifact type
    if artifact_type == ".png" or artifact_type == ".jpg" or artifact_type == ".jpeg":
        # Just copy the image file
        if output_dir:
            shutil.copy2(artifact_path, output_path)
        return artifact_path
        
    elif artifact_type == ".ply":
        # Render mesh or point cloud
        if "normal" in Path(artifact_path).name.lower():
            # Normal-estimated point cloud
            success = render_normal_cloud_to_image(artifact_path, output_path)
        else:
            # Regular mesh
            success = render_mesh_to_image(artifact_path, output_path)
            
        if success:
            return str(output_path)
        else:
            return None
            
    elif artifact_type == ".npy":
        # Assume it's a point cloud
        success = render_point_cloud_to_image(artifact_path, output_path)
        if success:
            return str(output_path)
        else:
            return None
            
    else:
        print(f"Unsupported artifact type: {artifact_type}")
        return None

def create_html_gallery(artifacts, output_dir, processed_images=None):
    """
    Create an HTML gallery of artifacts for browsing
    
    Args:
        artifacts: Dict of artifacts by structure and stage
        output_dir: Directory to store gallery files
        processed_images: Dict mapping artifact paths to image paths
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    if processed_images is None:
        processed_images = {}
    
    # Create HTML file
    html_path = output_dir / "gallery.html"
    
    # Begin HTML content
    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pipeline Artifacts Gallery</title>
    <style>
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, Cantarell, 'Open Sans', 'Helvetica Neue', sans-serif;
            line-height: 1.6;
            color: #333;
            margin: 0;
            padding: 20px;
            background-color: #f5f7fa;
        }
        
        h1, h2, h3 {
            color: #2c3e50;
        }
        
        h1 {
            border-bottom: 2px solid #3498db;
            padding-bottom: 10px;
            margin-bottom: 20px;
        }
        
        h2 {
            margin-top: 40px;
            border-bottom: 1px solid #ddd;
            padding-bottom: 5px;
        }
        
        .controls {
            margin-bottom: 20px;
            display: flex;
            gap: 15px;
            align-items: center;
            flex-wrap: wrap;
        }
        
        .filter-select {
            padding: 8px 12px;
            border: 1px solid #ddd;
            border-radius: 4px;
        }
        
        .gallery {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 40px;
        }
        
        .artifact-card {
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.05);
            transition: transform 0.2s;
        }
        
        .artifact-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 5px 15px rgba(0, 0, 0, 0.1);
        }
        
        .artifact-image {
            width: 100%;
            height: 200px;
            object-fit: contain;
            background: #f8f9fa;
            border-bottom: 1px solid #eee;
        }
        
        .artifact-info {
            padding: 15px;
        }
        
        .artifact-name {
            font-weight: bold;
            margin-bottom: 5px;
        }
        
        .artifact-path {
            font-size: 12px;
            color: #666;
            overflow: hidden;
            text-overflow: ellipsis;
            white-space: nowrap;
        }
        
        .artifact-stage {
            display: inline-block;
            background: #e1f5fe;
            color: #0277bd;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 12px;
            margin-top: 8px;
        }
        
        .structure-nav {
            position: sticky;
            top: 0;
            background: #fff;
            padding: 10px;
            border-bottom: 1px solid #ddd;
            margin-bottom: 20px;
            z-index: 100;
        }
        
        .structure-links {
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
        }
        
        .structure-link {
            text-decoration: none;
            padding: 5px 10px;
            background: #eee;
            color: #333;
            border-radius: 4px;
        }
        
        .structure-link:hover {
            background: #ddd;
        }
        
        .lightbox {
            display: none;
            position: fixed;
            z-index: 999;
            top: 0;
            left: 0;
            right: 0;
            bottom: 0;
            background: rgba(0,0,0,0.8);
            padding: 40px;
        }
        
        .lightbox-content {
            display: flex;
            flex-direction: column;
            align-items: center;
            justify-content: center;
            height: 100%;
        }
        
        .lightbox-image {
            max-width: 90%;
            max-height: 80%;
            object-fit: contain;
            background: white;
        }
        
        .lightbox-caption {
            color: white;
            margin-top: 20px;
            text-align: center;
        }
        
        .lightbox-close {
            position: absolute;
            top: 20px;
            right: 20px;
            color: white;
            font-size: 24px;
            cursor: pointer;
        }
        
        .lightbox-nav {
            position: absolute;
            top: 50%;
            transform: translateY(-50%);
            color: white;
            font-size: 30px;
            cursor: pointer;
            padding: 10px;
            background: rgba(0,0,0,0.5);
            border-radius: 50%;
            width: 30px;
            height: 30px;
            display: flex;
            align-items: center;
            justify-content: center;
        }
        
        .lightbox-prev {
            left: 20px;
        }
        
        .lightbox-next {
            right: 20px;
        }
    </style>
</head>
<body>
    <h1>Pipeline Artifacts Gallery</h1>
    
    <div class="controls">
        <select id="stageFilter" class="filter-select" onchange="filterArtifacts()">
            <option value="all">All Stages</option>
"""
    
    # Collect all stages
    all_stages = set()
    for structure_artifacts in artifacts.values():
        all_stages.update(structure_artifacts.keys())
    
    # Add stage options
    for stage in sorted(all_stages):
        html += f'            <option value="{stage}">{stage.replace("_", " ").title()}</option>\n'
    
    html += """        </select>
        
        <select id="typeFilter" class="filter-select" onchange="filterArtifacts()">
            <option value="all">All Artifact Types</option>
            <option value=".ply">Meshes (.ply)</option>
            <option value=".npy">Point Clouds (.npy)</option>
            <option value=".png">Images (.png)</option>
        </select>
    </div>
    
    <div class="structure-nav">
        <div class="structure-links">
"""
    
    # Add structure links
    for rcsb_id in sorted(artifacts.keys()):
        html += f'            <a href="#structure-{rcsb_id}" class="structure-link">{rcsb_id}</a>\n'
    
    html += """        </div>
    </div>
"""
    
    # Add structure sections
    for rcsb_id in sorted(artifacts.keys()):
        html += f'    <h2 id="structure-{rcsb_id}">{rcsb_id}</h2>\n'
        
        # Group artifacts by stage for this structure
        by_stage = artifacts[rcsb_id]
        
        for stage in sorted(by_stage.keys()):
            stage_artifacts = by_stage[stage]
            if not stage_artifacts:
                continue
                
            html += f'    <h3>{stage.replace("_", " ").title()}</h3>\n'
            html += '    <div class="gallery">\n'
            
            for artifact in stage_artifacts:
                artifact_path = artifact["path"]
                artifact_name = artifact["name"]
                artifact_type = artifact["type"]
                
                # Get image path for this artifact
                if artifact_path in processed_images:
                    image_path = processed_images[artifact_path]
                else:
                    # Use the path as is for images, otherwise indicate no preview
                    if artifact_type in [".png", ".jpg", ".jpeg"]:
                        image_path = artifact_path
                    else:
                        image_path = ""
                
                # Convert paths to relative if they're inside the output directory
                if image_path and str(output_dir) in image_path:
                    image_path = os.path.relpath(image_path, output_dir)
                
                # Only show artifacts with valid image representations
                if image_path:
                    html += f"""        <div class="artifact-card" data-stage="{stage}" data-type="{artifact_type}">
            <img 
                class="artifact-image" 
                src="{image_path}" 
                alt="{artifact_name}" 
                onclick="openLightbox(this)"
                data-path="{artifact_path}"
            >
            <div class="artifact-info">
                <div class="artifact-name">{artifact_name}</div>
                <div class="artifact-path">{artifact_path}</div>
                <span class="artifact-stage">{stage.replace("_", " ").title()}</span>
            </div>
        </div>
"""
            
            html += '    </div>\n'
    
    # Add lightbox for image viewing
    html += """    
    <!-- Lightbox for enlarged viewing -->
    <div class="lightbox" id="lightbox">
        <div class="lightbox-close" onclick="closeLightbox()">×</div>
        <div class="lightbox-nav lightbox-prev" onclick="navigateLightbox(-1)">←</div>
        <div class="lightbox-nav lightbox-next" onclick="navigateLightbox(1)">→</div>
        <div class="lightbox-content">
            <img class="lightbox-image" id="lightboxImage" src="">
            <div class="lightbox-caption" id="lightboxCaption"></div>
        </div>
    </div>
    
    <script>
        // Filtering functionality
        function filterArtifacts() {
            const stageFilter = document.getElementById('stageFilter').value;
            const typeFilter = document.getElementById('typeFilter').value;
            
            const cards = document.querySelectorAll('.artifact-card');
            
            cards.forEach(card => {
                const stage = card.getAttribute('data-stage');
                const type = card.getAttribute('data-type');
                
                let showCard = true;
                
                if (stageFilter !== 'all' && stage !== stageFilter) {
                    showCard = false;
                }
                
                if (typeFilter !== 'all' && !type.endsWith(typeFilter)) {
                    showCard = false;
                }
                
                card.style.display = showCard ? '' : 'none';
            });
        }
        
        // Lightbox functionality
        let currentIndex = 0;
        const visibleImages = [];
        
        function updateVisibleImages() {
            visibleImages.length = 0;
            document.querySelectorAll('.artifact-card:not([style*="display: none"]) .artifact-image').forEach(img => {
                visibleImages.push(img);
            });
        }
        
        function openLightbox(img) {
            updateVisibleImages();
            
            const lightbox = document.getElementById('lightbox');
            const lightboxImage = document.getElementById('lightboxImage');
            const lightboxCaption = document.getElementById('lightboxCaption');
            
            lightboxImage.src = img.src;
            lightboxCaption.textContent = img.getAttribute('data-path');
            
            // Find current index
            currentIndex = visibleImages.findIndex(image => image === img);
            
            lightbox.style.display = 'block';
        }
        
        function closeLightbox() {
            document.getElementById('lightbox').style.display = 'none';
        }
        
        function navigateLightbox(direction) {
            if (visibleImages.length === 0) return;
            
            currentIndex += direction;
            
            // Loop around if at the ends
            if (currentIndex < 0) currentIndex = visibleImages.length - 1;
            if (currentIndex >= visibleImages.length) currentIndex = 0;
            
            const img = visibleImages[currentIndex];
            const lightboxImage = document.getElementById('lightboxImage');
            const lightboxCaption = document.getElementById('lightboxCaption');
            
            lightboxImage.src = img.src;
            lightboxCaption.textContent = img.getAttribute('data-path');
        }
        
        // Key navigation for lightbox
        document.addEventListener('keydown', function(e) {
            if (document.getElementById('lightbox').style.display === 'block') {
                if (e.key === 'Escape') {
                    closeLightbox();
                } else if (e.key === 'ArrowLeft') {
                    navigateLightbox(-1);
                } else if (e.key === 'ArrowRight') {
                    navigateLightbox(1);
                }
            }
        });
    </script>
</body>
</html>
"""
    
    # Write HTML file
    with open(html_path, 'w') as f:
        f.write(html)
    
    print(f"Gallery created at {html_path}")
    return html_path

def create_artifact_viewer(artifacts, output_dir=None, max_structures=20, gallery=True, parallel=False):
    """
    Create a viewer for the selected artifacts
    
    Args:
        artifacts: Dict of artifacts by structure and stage
        output_dir: Directory to store generated visualizations
        max_structures: Maximum number of structures to process
        gallery: Whether to create an HTML gallery
        parallel: Whether to process artifacts in parallel
    """
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
    
    # Limit number of structures
    if len(artifacts) > max_structures:
        print(f"Limiting to {max_structures} structures (from {len(artifacts)})")
        struct_ids = list(artifacts.keys())[:max_structures]
        limited_artifacts = {id: artifacts[id] for id in struct_ids}
        artifacts = limited_artifacts
    
    # Flatten artifacts list for processing
    flat_artifacts = []
    for rcsb_id, stages in artifacts.items():
        for stage, stage_artifacts in stages.items():
            flat_artifacts.extend(stage_artifacts)
    
    print(f"Processing {len(flat_artifacts)} artifacts for {len(artifacts)} structures")
    
    # Create output directory if needed
    if gallery and not output_dir:
        output_dir = Path("pipeline_artifacts_gallery")
        output_dir.mkdir(exist_ok=True, parents=True)
        print(f"Created gallery directory: {output_dir}")
    
    # Process artifacts
    processed_images = {}
    
    if parallel and len(flat_artifacts) > 5:
        print(f"Processing artifacts in parallel...")
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            # Submit all tasks
            future_to_artifact = {
                executor.submit(process_artifact, artifact, output_dir): artifact
                for artifact in flat_artifacts 
                if artifact['exists']
            }
            
            # Process results as they complete
            for i, future in enumerate(concurrent.futures.as_completed(future_to_artifact)):
                artifact = future_to_artifact[future]
                try:
                    image_path = future.result()
                    if image_path:
                        processed_images[artifact['path']] = image_path
                except Exception as e:
                    print(f"Error processing {artifact['name']}: {e}")
                
                # Show progress
                print(f"Processed {i+1}/{len(future_to_artifact)} artifacts", end="\r")
            print()  # New line after progress
    
    # Create gallery if requested
    if gallery and output_dir:
        gallery_path = create_html_gallery(artifacts, output_dir, processed_images)
        print(f"Opening gallery in browser: {gallery_path}")
        webbrowser.open(f"file://{os.path.abspath(gallery_path)}")
        return gallery_path
    
    return processed_images

def main():
    """Main function"""
    args = parse_args()
    
    # Get structures to process
    structures = []
    
    if args.structures:
        # Use structures provided on command line
        structures = args.structures
    else:
        # Prompt user to enter or paste structures
        print("Enter structure IDs (one per line, blank line to finish):")
        while True:
            line = input()
            if not line:
                break
            # Split in case multiple space-separated IDs were pasted
            structures.extend(line.split())
    
    if not structures:
        print("No structures specified. Exiting.")
        return
    
    # Find logs for structures
    structure_logs = find_structure_logs(structures, args.log_dir)
    
    if not structure_logs:
        print("No logs found for specified structures.")
        return
    
    print(f"Found logs for {len(structure_logs)} out of {len(structures)} structures")
    
    # Collect artifacts
    artifacts = collect_artifacts(structure_logs, args.stage)
    
    # Filter artifacts for visualization
    vis_artifacts = filter_artifacts_by_type(artifacts, file_type=args.type)
    
    # Create viewer for artifacts
    output_dir = args.output_dir if args.output_dir else "pipeline_gallery"
    
    create_artifact_viewer(
        vis_artifacts, 
        output_dir, 
        args.max_structures,
        args.gallery,
        args.parallel
    )

if __name__ == "__main__":
    main()