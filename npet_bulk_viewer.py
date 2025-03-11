import os
import sys
import glob
import json
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QLineEdit, QLabel, QGridLayout, QScrollArea, QCheckBox
from PyQt5.QtCore import Qt
import Bio.PDB
from Bio.PDB.PDBParser import PDBParser
from pathlib import Path
from scipy.spatial.distance import cdist

from ribctl.lib.schema.types_ribosome import ConstrictionSite, PTCInfo
from ribctl.asset_manager.asset_types import AssetType
from ribctl.ribosome_ops import RibosomeOps


class EnhancedMeshViewer(QWidget):
    def __init__(self):
        super().__init__()
        self.meshes = {}
        self.plotter = None
        self.grid_plotters = []
        self.current_grid_page = 0
        self.meshes_per_page = 9  # 3x3 grid
        self.show_landmarks = True
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle('Enhanced PLY Mesh Viewer')
        self.setGeometry(100, 100, 1200, 800)
        
        main_layout = QVBoxLayout()
        
        # Control panel
        control_panel = QHBoxLayout()
        
        # Input for mesh directory
        self.dir_label = QLabel('Mesh Directory:')
        self.dir_input = QLineEdit()
        self.dir_input.setText(os.getcwd())  # Default to current directory
        self.load_btn = QPushButton('Load Meshes')
        self.load_btn.clicked.connect(self.load_meshes)
        
        # Input for specific mesh ID
        self.id_label = QLabel('Mesh ID:')
        self.id_input = QLineEdit()
        self.view_btn = QPushButton('View Single Mesh')
        self.view_btn.clicked.connect(self.view_single_mesh)
        
        # Grid navigation
        self.prev_btn = QPushButton('Previous Page')
        self.prev_btn.clicked.connect(self.prev_page)
        self.next_btn = QPushButton('Next Page')
        self.next_btn.clicked.connect(self.next_page)
        self.page_label = QLabel('Page: 1')
        
        # Toggle landmarks
        self.landmarks_check = QCheckBox('Show Landmarks')
        self.landmarks_check.setChecked(True)
        self.landmarks_check.stateChanged.connect(self.toggle_landmarks)
        
        # Add widgets to control panel
        control_panel.addWidget(self.dir_label)
        control_panel.addWidget(self.dir_input)
        control_panel.addWidget(self.load_btn)
        control_panel.addWidget(self.id_label)
        control_panel.addWidget(self.id_input)
        control_panel.addWidget(self.view_btn)
        control_panel.addWidget(self.prev_btn)
        control_panel.addWidget(self.page_label)
        control_panel.addWidget(self.next_btn)
        control_panel.addWidget(self.landmarks_check)
        
        main_layout.addLayout(control_panel)
        
        # Grid layout for mesh thumbnails
        self.grid_container = QWidget()
        self.grid_layout = QGridLayout(self.grid_container)
        
        # Scroll area for grid
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(self.grid_container)
        
        main_layout.addWidget(scroll_area)
        
        self.setLayout(main_layout)
        
    def toggle_landmarks(self, state):
        self.show_landmarks = (state == Qt.Checked)
        self.update_grid()
        
    def load_meshes(self):
        # Clear existing data
        self.meshes = {}
        for plotter in self.grid_plotters:
            plotter.close()
        self.grid_plotters = []
        
        # Clear grid layout
        for i in reversed(range(self.grid_layout.count())):
            widget = self.grid_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)
        
        # Get mesh directory
        mesh_dir = self.dir_input.text()
        mesh_files = glob.glob(os.path.join(mesh_dir, '*_NPET_MESH.ply'))
        
        if not mesh_files:
            print("No PLY mesh files found in the specified directory.")
            return
            
        # Load meshes
        for file_path in mesh_files:
            mesh_id = os.path.basename(file_path).split('_')[0]
            self.meshes[mesh_id] = file_path
        
        print(f"Loaded {len(self.meshes)} meshes.")
        self.current_grid_page = 0
        self.update_grid()
    
    def update_grid(self):
        # Close existing plotters
        for plotter in self.grid_plotters:
            plotter.close()
        self.grid_plotters = []
        
        # Clear grid layout
        for i in reversed(range(self.grid_layout.count())): 
            widget = self.grid_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)
        
        # Get mesh IDs for the current page
        mesh_ids = list(self.meshes.keys())
        start_idx = self.current_grid_page * self.meshes_per_page
        end_idx = min(start_idx + self.meshes_per_page, len(mesh_ids))
        current_mesh_ids = mesh_ids[start_idx:end_idx]
        
        # Calculate grid dimensions for a 3x3 grid
        grid_size = 3
        
        # Create grid of mesh viewers
        for idx, mesh_id in enumerate(current_mesh_ids):
            row = (idx // grid_size) * 2  # Multiply by 2 to leave room for labels
            col = idx % grid_size
            
            # Create a small plotter for each mesh
            plotter = QtInteractor()
            plotter.setMinimumSize(300, 300)
            
            # Load and display mesh
            mesh_file = self.meshes[mesh_id]
            mesh = pv.read(mesh_file)
            plotter.add_mesh(mesh, color='lightblue', show_edges=True, opacity=0.7)
            
            # Add landmarks if enabled
            if self.show_landmarks:
                try:
                    self.add_landmarks_to_plotter(plotter, mesh_id)
                except Exception as e:
                    print(f"Error adding landmarks for {mesh_id}: {e}")
                    
            plotter.add_text(mesh_id, position='upper_left', font_size=10)
            plotter.reset_camera()
            
            # Add plotter to grid
            self.grid_layout.addWidget(plotter, row, col)
            self.grid_plotters.append(plotter)
            
            # Add a label with the mesh ID
            id_label = QLabel(f"ID: {mesh_id}")
            id_label.setAlignment(Qt.AlignCenter)
            id_label.setStyleSheet("background-color: rgba(200, 200, 200, 150); padding: 5px;")
            
            # Add the label to the grid layout in the cell below the plotter
            self.grid_layout.addWidget(id_label, row + 1, col)
        
        # Update page label
        total_pages = max(1, int(np.ceil(len(self.meshes) / self.meshes_per_page)))
        self.page_label.setText(f'Page: {self.current_grid_page + 1}/{total_pages}')
        
        # Enable/disable navigation buttons
        self.prev_btn.setEnabled(self.current_grid_page > 0)
        self.next_btn.setEnabled(self.current_grid_page < total_pages - 1)
    
    def add_landmarks_to_plotter(self, plotter, mesh_id):
        """Add PTC, constriction site, and protein segments to the plotter"""
        try:
            # Initialize RibosomeOps for access to structure data
            ro = RibosomeOps(mesh_id)
            
            # Load PTC information
            ptc_path = AssetType.PTC.get_path(mesh_id)
            if ptc_path.exists():
                with open(ptc_path, 'r') as f:
                    ptc_info = PTCInfo.model_validate(json.load(f))
                ptc_location = np.array(ptc_info.location)
                # Add PTC as a yellow sphere
                plotter.add_mesh(pv.Sphere(center=ptc_location, radius=5), color='yellow')
            
            # Load constriction site information
            constriction_path = AssetType.CONSTRICTION_SITE.get_path(mesh_id)
            if constriction_path.exists():
                with open(constriction_path, 'r') as f:
                    constriction_info = ConstrictionSite.model_validate(json.load(f))
                constriction_location = np.array(constriction_info.location)
                # Add constriction site as a red sphere
                plotter.add_mesh(pv.Sphere(center=constriction_location, radius=5), color='red')
            else:
                # If no constriction site file exists, we can't find the closest residues
                return
            
            # Find uL4 and uL22 protein chains
            profile = ro.profile
            ul4_chain = None
            ul22_chain = None
            
            # Search through proteins to find uL4 and uL22
            for protein in profile.proteins:
                for nomenclature in protein.nomenclature:
                    if nomenclature.value == "uL4":
                        ul4_chain = protein.auth_asym_id
                    elif nomenclature.value == "uL22":
                        ul22_chain = protein.auth_asym_id
            
            # If we found the chains, get the residues closest to constriction site
            if ul4_chain or ul22_chain:
                model = ro.assets.biopython_structure()[0]
                
                # Process uL4 chain
                if ul4_chain and ul4_chain in model.child_dict:
                    ul4_residues = list(model[ul4_chain].get_residues())
                    ul4_points = self.get_closest_residues(ul4_residues, constriction_location, 15)
                    if ul4_points.size > 0:
                        # Add uL4 residues as blue points
                        plotter.add_mesh(pv.PolyData(ul4_points), color='blue', point_size=10, render_points_as_spheres=True)
                
                # Process uL22 chain
                if ul22_chain and ul22_chain in model.child_dict:
                    ul22_residues = list(model[ul22_chain].get_residues())
                    ul22_points = self.get_closest_residues(ul22_residues, constriction_location, 15)
                    if ul22_points.size > 0:
                        # Add uL22 residues as orange points
                        plotter.add_mesh(pv.PolyData(ul22_points), color='orange', point_size=10, render_points_as_spheres=True)
        
        except Exception as e:
            print(f"Error processing landmarks for {mesh_id}: {str(e)}")
        
    def get_closest_residues(self, residues, target_location, n=15):
        """Get the n residues with centers closest to the target location"""
        residue_centers = []
        valid_residues = []
        
        for residue in residues:
            try:
                # Calculate center of residue
                coords = np.array([atom.coord for atom in residue.get_atoms()])
                if len(coords) > 0:
                    center = coords.mean(axis=0)
                    residue_centers.append(center)
                    valid_residues.append(residue)
            except Exception:
                continue
        
        if not residue_centers:
            return np.array([])
        
        # Calculate distances to target location
        residue_centers = np.array(residue_centers)
        distances = cdist([target_location], residue_centers)[0]
        
        # Get indices of n closest residues
        closest_indices = np.argsort(distances)[:n]
        
        # Convert residues to points (use CA atom or average position)
        points = []
        for idx in closest_indices:
            try:
                residue = valid_residues[idx]
                if 'CA' in residue:
                    points.append(residue['CA'].coord)
                else:
                    coords = np.array([atom.coord for atom in residue.get_atoms()])
                    points.append(coords.mean(axis=0))
            except Exception:
                continue
        
        return np.array(points)
    
    def prev_page(self):
        if self.current_grid_page > 0:
            self.current_grid_page -= 1
            self.update_grid()
    
    def next_page(self):
        total_pages = int(np.ceil(len(self.meshes) / self.meshes_per_page))
        if self.current_grid_page < total_pages - 1:
            self.current_grid_page += 1
            self.update_grid()
    
    def open_mesh_by_id(self, mesh_id):
        """Helper method to open a mesh by its ID"""
        self.id_input.setText(mesh_id)
        self.view_single_mesh()
    
    def view_single_mesh(self):
        mesh_id = self.id_input.text().strip()
        
        if not mesh_id:
            print("Please enter a mesh ID.")
            return
        
        if mesh_id not in self.meshes:
            print(f"Mesh ID '{mesh_id}' not found.")
            return
        
        # Close previous plotter if exists
        if self.plotter is not None:
            self.plotter.close()
        
        # Create new plotter window
        self.plotter = pv.Plotter()
        
        # Load and display mesh
        mesh = pv.read(self.meshes[mesh_id])
        self.plotter.add_mesh(mesh, color='lightblue', show_edges=True, opacity=0.7)
        
        # Add landmarks
        if self.show_landmarks:
            try:
                self.add_landmarks_to_plotter(self.plotter, mesh_id)
            except Exception as e:
                print(f"Error adding landmarks: {e}")
        
        # Add mesh information
        info_text = [
            f"Mesh ID: {mesh_id}",
            f"Vertices: {mesh.n_points}",
            f"Faces: {mesh.n_faces}",
            f"File: {os.path.basename(self.meshes[mesh_id])}"
        ]
        self.plotter.add_text('\n'.join(info_text), position='upper_left', font_size=12)
        
        # Show the plotter
        self.plotter.show(title=f"Mesh Viewer - {mesh_id}")


if __name__ == '__main__':
    try:
        app = QApplication(sys.argv)
        viewer = EnhancedMeshViewer()
        viewer.show()
        sys.exit(app.exec_())
    except ModuleNotFoundError as e:
        if 'pyvistaqt' in str(e):
            print("Error: Missing pyvistaqt package. Please install it with:")
            print("pip install pyvistaqt")
        else:
            print(f"Error: {e}")
            print("You might need to install the required packages:")
            print("pip install pyvista pyvistaqt numpy PyQt5 scipy biopython")