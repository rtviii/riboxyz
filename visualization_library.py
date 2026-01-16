"""
NPET Visualization Function Registry
"""
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from typing import Optional

# -------------------------------------------------------------------
## --- REGISTRABLE VISUALIZATION FUNCTIONS ---
#
# These are the refactored functions. They take a `plotter`
# and *already loaded* data. They just add actors.
# -------------------------------------------------------------------

def viz_simple_mesh(
    plotter: QtInteractor, 
    mesh_to_show: pv.PolyData, 
    rcsb_id: str
):
    """
    Adds a simple mesh and a text label to an existing plotter.
    """
    if mesh_to_show is None:
        plotter.add_text("Artifact file not found.", position='upper_left', color='red')
        return

    plotter.add_mesh(
        mesh_to_show, 
        color='lightblue', 
        show_edges=True, 
        opacity=0.8
    )
    plotter.add_text(
        f'RCSB_ID: {rcsb_id}',
        position='upper_right', 
        font_size=10, 
        font='courier'
    )
    plotter.reset_camera()

def viz_single_landmark(
    plotter: QtInteractor,
    landmark_json: dict,
    rcsb_id: str
):
    """
    Visualizes a single landmark from its JSON file.
    """
    if landmark_json is None:
        plotter.add_text("Landmark JSON not found.", position='upper_left', color='red')
        return
        
    try:
        coord = np.array(landmark_json['location'])
        plotter.add_points(
            coord,
            color='red',
            point_size=20,
            render_points_as_spheres=True,
            label="Landmark"
        )
    except KeyError:
        plotter.add_text("Landmark JSON missing 'location' key.", position='upper_left', color='red')
    except Exception as e:
        plotter.add_text(f"Error plotting landmark: {e}", position='upper_left', color='red')

    plotter.add_text(
        f'RCSB_ID: {rcsb_id}',
        position='upper_right', 
        font_size=10, 
        font='courier'
    )
    plotter.reset_camera()

def viz_landmarks_and_pointcloud(
    plotter: QtInteractor,
    rcsb_id: str,
    background_mesh: Optional[pv.PolyData] = None,
    point_cloud_data: Optional[np.ndarray] = None,
    ptc_coord: Optional[np.ndarray] = None,
    constriction_coord: Optional[np.ndarray] = None,
):
    """
    Visualizes landmarks and/or a point cloud,
    optionally on top of a background mesh.
    """
    
    # 1. Add the background mesh if provided
    if background_mesh is not None:
        plotter.add_mesh(
            background_mesh,
            style='surface',
            color='white',
            opacity=0.1,
            show_edges=False
        )

    # 2. Add the main point cloud if provided
    if point_cloud_data is not None:
        cloud = pv.PolyData(point_cloud_data)
        plotter.add_mesh(
            cloud, 
            color='grey', 
            point_size=3, 
            opacity=0.1,
            style='points_gaussian',
            emissive=True
        )
    
    # 3. Add PTC landmark if provided
    if ptc_coord is not None:
        plotter.add_points(
            ptc_coord,
            color='red',
            point_size=15,
            render_points_as_spheres=True,
            label="PTC"
        )
        
    # 4. Add Constriction landmark if provided
    if constriction_coord is not None:
        plotter.add_points(
            constriction_coord,
            color='blue',
            point_size=15,
            render_points_as_spheres=True,
            label="Constriction"
        )

    plotter.add_legend(bcolor=None)
    plotter.add_text(
        f'RCSB_ID: {rcsb_id}',
        position='upper_right', 
        font_size=10, 
        font='courier'
    )
    plotter.reset_camera()

import gemmi

def viz_ribosome_structure(
    plotter: QtInteractor,
    cif_path: str,
    rcsb_id: str,
    chains_to_show: Optional[list] = None,
    alpha: float = 0.5,
    point_size: int = 5
):
    """
    Visualizes ribosome structure as residue centers from CIF file.
    Uses BioPython for fast parsing.
    """
    try:
        from Bio.PDB import MMCIFParser
        import warnings
        warnings.filterwarnings('ignore')  # Suppress BioPython warnings
        
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure(rcsb_id, cif_path)
        
        # Get first model
        model = structure[0]
        
        # Collect residue centers of mass quickly
        residue_centers = []
        for chain in model:
            if chains_to_show and chain.id not in chains_to_show:
                continue
            for residue in chain:
                # Get all atoms in residue
                atoms = [atom for atom in residue.get_atoms()]
                if atoms:
                    # Quick center of mass calculation
                    coords = np.array([atom.coord for atom in atoms])
                    center = coords.mean(axis=0)
                    residue_centers.append(center)
        
        if not residue_centers:
            plotter.add_text("No residues found", position='lower_left', color='yellow')
            return
            
        residue_centers = np.array(residue_centers)
        
        # Add as point cloud with better visibility
        cloud = pv.PolyData(residue_centers)
        plotter.add_mesh(
            cloud,
            name='ribosome_structure',
            color='lime',  # Brighter green
            point_size=point_size,
            opacity=alpha,
            render_points_as_spheres=True,
            label="Ribosome Structure",
            emissive=True  # Makes points glow slightly
        )
        
        plotter.add_text(
            f"Ribosome: {len(residue_centers)} residues",
            position='lower_left',
            font_size=10,
            color='lime',
            name='ribosome_label'
        )
        
    except Exception as e:
        plotter.add_text(f"Error loading structure: {e}", position='lower_left', color='red')
        import traceback
        print(traceback.format_exc())