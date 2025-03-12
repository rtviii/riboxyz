#!/usr/bin/env python
import os
import sys
import json
import argparse
import numpy as np
import trimesh
import pyrender
import matplotlib.pyplot as plt
from pathlib import Path


def load_ply_mesh(mesh_path):
    """Load a PLY mesh file using trimesh"""
    try:
        return trimesh.load(mesh_path)
    except Exception as e:
        print(f"Error loading mesh from {mesh_path}: {e}")
        return None


def load_landmarks(landmarks_path):
    """Load landmarks from JSON file"""
    try:
        with open(landmarks_path, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Error loading landmarks from {landmarks_path}: {e}")
        return None


def create_sphere(center, radius, color):
    """Create a sphere mesh at the specified center point"""
    sphere = trimesh.creation.icosphere(radius=radius)
    sphere.apply_translation(center)
    sphere.visual.face_colors = color
    return sphere


def create_point_cloud(points, color):
    """Create a point cloud from a list of points"""
    if not points or len(points) == 0:
        return None
    
    points = np.array(points)
    cloud = trimesh.points.PointCloud(points)
    cloud.colors = np.tile(color, (len(points), 1))
    return cloud


def visualize_structure(rcsb_id, mesh_dir, landmarks_dir, output_dir=None):
    """Visualize a structure with its landmarks using Trimesh"""
    # Build file paths
    mesh_path = os.path.join(mesh_dir, f"{rcsb_id}_NPET_MESH.ply")
    landmarks_path = os.path.join(landmarks_dir, f"{rcsb_id}_landmarks_pts.json")
    
    # Check if files exist
    if not os.path.exists(mesh_path):
        print(f"Mesh file not found: {mesh_path}")
        return False
    
    if not os.path.exists(landmarks_path):
        print(f"Landmarks file not found: {landmarks_path}")
        return False
    
    # Load mesh and landmarks
    mesh = load_ply_mesh(mesh_path)
    landmarks = load_landmarks(landmarks_path)
    
    if mesh is None or landmarks is None:
        return False
    
    # Create a scene
    scene = trimesh.Scene()
    
    # Add the tunnel mesh
    mesh.visual.face_colors = [100, 100, 100, 150]  # Gray with some transparency
    scene.add_geometry(mesh, node_name='tunnel')
    
    # Add PTC landmark (yellow sphere)
    if landmarks.get("ptc"):
        ptc_sphere = create_sphere(landmarks["ptc"], radius=5.0, color=[255, 255, 0, 255])
        scene.add_geometry(ptc_sphere, node_name='ptc')
    
    # Add constriction site landmark (red sphere)
    if landmarks.get("constriction"):
        constriction_sphere = create_sphere(landmarks["constriction"], radius=5.0, color=[255, 0, 0, 255])
        scene.add_geometry(constriction_sphere, node_name='constriction')
    
    # Add uL4 protein residues (blue points)
    if landmarks.get("uL4") and len(landmarks["uL4"]) > 0:
        ul4_cloud = create_point_cloud(landmarks["uL4"], color=[0, 0, 255, 255])
        if ul4_cloud:
            scene.add_geometry(ul4_cloud, node_name='uL4')
    
    # Add uL22 protein residues (orange points)
    if landmarks.get("uL22") and len(landmarks["uL22"]) > 0:
        ul22_cloud = create_point_cloud(landmarks["uL22"], color=[255, 165, 0, 255])
        if ul22_cloud:
            scene.add_geometry(ul22_cloud, node_name='uL22')
    
    # Create output directory if specified
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"{rcsb_id}_visualization.png")
        
        # Save the visualization to an image file
        try:
            # Get a good camera angle
            camera_params = scene.camera_transform
            scene.camera_transform = trimesh.transformations.rotation_matrix(
                angle=np.radians(30),
                direction=[0, 1, 0],
                point=scene.centroid
            )
            
            # Render the scene to an image
            png = scene.save_image(resolution=[1920, 1080], visible=True)
            with open(output_path, 'wb') as f:
                f.write(png)
            print(f"Saved visualization to {output_path}")
            
            # Also create an interactive visualization
            interactive_output = os.path.join(output_dir, f"{rcsb_id}_interactive.html")
            scene.export(interactive_output)
            print(f"Saved interactive visualization to {interactive_output}")
        except Exception as e:
            print(f"Error saving visualization: {e}")
    
    # Show the visualization interactively
    scene.show()
    
    return True


def main():
    parser = argparse.ArgumentParser(description='Visualize a ribosome structure with landmarks using Trimesh')
    parser.add_argument('rcsb_id', help='RCSB PDB ID of the structure to visualize')
    parser.add_argument('--mesh-dir', required=True, help='Directory containing NPET mesh files')
    parser.add_argument('--landmarks-dir', required=True, help='Directory containing landmark JSON files')
    parser.add_argument('--output-dir', help='Directory to save visualization images (optional)')
    args = parser.parse_args()
    
    # Visualize the structure
    success = visualize_structure(
        args.rcsb_id, 
        args.mesh_dir, 
        args.landmarks_dir, 
        args.output_dir
    )
    
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()