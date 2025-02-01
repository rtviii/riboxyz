# import warnings
# from Bio.PDB.MMCIFParser import MMCIFParser
# import numpy as np
# import random
# import os
# import alphashape
# import trimesh
# from ribctl.asset_manager.asset_types import AssetType

# data_dir = os.environ.get("DATA_DIR")
# warnings.filterwarnings("ignore")


# def clean_mesh(mesh):
#     """Clean up a mesh by removing duplicate vertices and faces."""
#     if mesh is None:
#         print("Warning: Received None mesh in clean_mesh")
#         return None

#     try:
#         # Remove duplicate vertices
#         mesh.merge_vertices()

#         # Remove duplicate faces
#         mesh.remove_duplicate_faces()

#         # Remove degenerate faces
#         mesh.remove_degenerate_faces()

#         # Ensure consistent face winding
#         mesh.fix_normals()

#         return mesh
#     except Exception as e:
#         print(f"Error in clean_mesh: {e}")
#         return mesh

# def sample_within_alpha_shape(alpha_shape_mesh, num_samples):
#     """Samples points within the boundaries of a 3D alpha shape."""
#     if alpha_shape_mesh is None:
#         raise ValueError("Received None mesh in sample_within_alpha_shape")

#     bbox_min, bbox_max = alpha_shape_mesh.bounds

#     sampled_points = []
#     max_attempts = 10  # Prevent infinite loops
#     attempt = 0

#     while len(sampled_points) < num_samples and attempt < max_attempts:
#         points = np.random.uniform(
#             low=bbox_min, high=bbox_max, size=(num_samples * 4, 3)
#         )
#         inside = alpha_shape_mesh.contains(points)
#         sampled_points.extend(points[inside])
#         print(f"Attempt {attempt + 1}: Found {len(sampled_points)} points")
#         attempt += 1

#     if len(sampled_points) == 0:
#         raise ValueError("No points sampled within alpha shape - mesh may be invalid")

#     sampled_points = np.array(sampled_points[:num_samples])
#     print(f"Final number of sampled points: {len(sampled_points)}")
#     return sampled_points


# def save_alpha_shape_as_ply(alpha_shape, file_path):
#     """Save a watertight Trimesh alpha shape as a PLY file."""
#     if alpha_shape is None:
#         raise ValueError("Cannot save None mesh")

#     if not isinstance(alpha_shape, trimesh.Trimesh):
#         raise TypeError("alpha_shape must be a trimesh.Trimesh object")

#     # Check mesh orientation and flip if needed
#     if alpha_shape.volume < 0:
#         print("Negative volume detected, flipping mesh orientation...")
#         alpha_shape.invert()

#     # Add validation
#     if not alpha_shape.is_watertight:
#         print("Mesh is not watertight, attempting to fix...")
#         try:
#             # Try to fix holes while preserving the mesh object
#             vertices_to_fix = alpha_shape.fill_holes()
#             if vertices_to_fix is not None and isinstance(vertices_to_fix, bool):
#                 print("Warning: Could not fix holes in mesh")
#             # Additional cleanup attempts
#             alpha_shape.remove_degenerate_faces()
#             alpha_shape.remove_duplicate_faces()
#             alpha_shape.remove_unreferenced_vertices()
#         except Exception as e:
#             print(f"Warning: Error while trying to fix holes: {e}")

#     # Check winding and fix if needed
#     try:
#         if not alpha_shape.is_winding_consistent:
#             print("Fixing inconsistent winding...")
#             alpha_shape.fix_normals()
#     except Exception as e:
#         print(f"Warning: Error while checking winding consistency: {e}")

#     try:
#         alpha_shape.export(file_path, file_type="ply", encoding="ascii")
#         print(f"Alpha shape saved as PLY file at {file_path}")
#     except Exception as e:
#         print(f"Error saving PLY file: {e}")
#         raise

#     return alpha_shape


# def validate_mesh_pyvista(mesh, stage="unknown"):
#     """Validate and print mesh properties."""
#     if mesh is None:
#         print(f"WARNING: Null mesh at stage {stage}")
#         return None

#     print(f"\nMesh properties at stage: {stage}")
#     print(f"- Is watertight: {mesh.is_watertight}")
#     print(f"- Is winding consistent: {mesh.is_winding_consistent}")
#     print(f"- Number of vertices: {len(mesh.vertices)}")
#     print(f"- Number of faces: {len(mesh.faces)}")
#     print(f"- Volume: {mesh.volume}")
#     print(f"- Surface area: {mesh.area}")
#     return mesh


# def produce_alpha_contour(RCSB_ID, alpha):
#     """Produce a watertight alpha shape contour with normal orientation check."""

#     num_samples = 5000
#     cifpath = AssetType.MMCIF.get_path(RCSB_ID)

#     # Get initial point cloud
#     point_cloud = cif_to_point_cloud(cifpath)
#     point_cloud = np.array(point_cloud)
#     print(f"Initial point cloud size: {len(point_cloud)}")

#     if len(point_cloud) == 0:
#         raise ValueError(f"Empty point cloud for {RCSB_ID}")

#     # Create initial alpha shape
#     print(f"Constructing alpha shape with a={alpha}")
#     initial_shape = alphashape.alphashape(point_cloud, alpha)
#     validate_mesh_pyvista(initial_shape, "initial_shape")

#     # Split into components
#     components = initial_shape.split(only_watertight=False)
#     if not components:
#         raise ValueError("No components found in alpha shape")

#     print(f"Found {len(components)} components")

#     # Get largest component
#     alpha_shape = max(components, key=lambda c: abs(c.volume))
#     alpha_shape = clean_mesh(alpha_shape)
#     validate_mesh_pyvista(alpha_shape, "after_clean")

#     # Check orientation and flip if needed
#     if alpha_shape.volume < 0:
#         print("Flipping mesh orientation...")
#         alpha_shape.invert()
#         validate_mesh_pyvista(alpha_shape, "after_flip")

#     # Sample points
#     random.seed(10)
#     print(f"\nResampling {num_samples} points")
#     new_points = sample_within_alpha_shape(alpha_shape, num_samples)

#     if len(new_points) == 0:
#         raise ValueError("No points sampled from alpha shape")

#     # Create final alpha shape
#     alpha_shape_renew = alphashape.alphashape(new_points, alpha)
#     alpha_shape_renew = clean_mesh(alpha_shape_renew)
#     validate_mesh_pyvista(alpha_shape_renew, "final_shape")

#     # Final orientation check
#     if alpha_shape_renew.volume < 0:
#         print("Flipping final mesh orientation...")
#         alpha_shape_renew.invert()
#         validate_mesh_pyvista(alpha_shape_renew, "after_final_flip")

#     path = AssetType.ALPHA_SHAPE.get_path(RCSB_ID)
#     alpha_shape_renew = save_alpha_shape_as_ply(alpha_shape_renew, path)
#     print(f"Saved to {path}")

#     return alpha_shape_renew
