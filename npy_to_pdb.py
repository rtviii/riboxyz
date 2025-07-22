def npy_to_pdb(input_npy_path: str, output_pdb_path: str, alpha_shape_path: str = None) -> str:
    """
    Convert a numpy array of 3D coordinates to a PDB file, optionally clipping with alpha shape.
    
    Args:
        input_npy_path: Path to .npy file containing Nx3 array of coordinates
        output_pdb_path: Path where PDB file will be written
        alpha_shape_path: Optional path to alpha shape mesh file for clipping
        
    Returns:
        str: Path to the generated PDB file
    """
    import numpy as np
    import pyvista as pv
    from pathlib import Path
    
    # Load the coordinates
    coords = np.load(input_npy_path)
    
    # Validate input
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"Expected Nx3 array, got shape {coords.shape}")
    
    print(f"Loaded {len(coords)} points from {input_npy_path}")
    
    # Clip with alpha shape if provided
    if alpha_shape_path:
        print(f"Clipping with alpha shape: {alpha_shape_path}")
        
        # Load alpha shape mesh
        alpha_mesh = pv.read(alpha_shape_path)
        
        # Create point cloud and select points inside mesh
        point_cloud = pv.PolyData(coords)
        selection = point_cloud.select_enclosed_points(alpha_mesh, check_surface=False)
        mask = selection["SelectedPoints"]
        
        # Keep only interior points
        coords = coords[mask == 1]
        print(f"After clipping: {len(coords)} points remain")
        
        if len(coords) == 0:
            raise ValueError("No points remain after alpha shape clipping!")
    
    # Get filename for header
    filename = Path(input_npy_path).stem
    
    with open(output_pdb_path, 'w') as pdb_file:
        # Write PDB header
        pdb_file.write(f"HEADER    COORDINATES FROM {filename.upper()}\n")
        pdb_file.write(f"TITLE     CONVERTED FROM {input_npy_path}\n")
        if alpha_shape_path:
            pdb_file.write(f"REMARK    CLIPPED WITH {alpha_shape_path}\n")
        pdb_file.write("REMARK    Generated from numpy array\n")
        
        atom_num = 1
        connect_records = []
        
        # Write each coordinate as a CA atom in chain A, residue 1
        for point_idx, (x, y, z) in enumerate(coords):
            pdb_file.write(
                f"ATOM  {atom_num:5d}  CA  ALA A   1    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n"
            )
            
            # Create bonds between consecutive atoms
            if point_idx > 0:
                connect_records.append(f"CONECT{atom_num-1:5d}{atom_num:5d}\n")
            
            atom_num += 1
        
        # Write chain terminator
        pdb_file.write(f"TER   {atom_num:5d}      ALA A   1\n")
        
        # Write all CONECT records
        for connect in connect_records:
            pdb_file.write(connect)
            
        # Write end record
        pdb_file.write("END\n")
    
    print(f"Generated PDB: {output_pdb_path}")
    return output_pdb_path




npy_to_pdb("/Users/rtviii/dev/RIBETL_DATA/4UG0/artifacts/4UG0_refined_cluster.npy", "/Users/rtviii/dev/RIBETL_DATA/4UG0/artifacts/4UG0_refined_cluster.pdb",
            '/Users/rtviii/dev/RIBETL_DATA/4UG0/4UG0_NPET_MESH_ascii.ply')