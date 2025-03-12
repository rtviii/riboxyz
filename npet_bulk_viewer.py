#!/usr/bin/env python
import os
import sys
import glob
import json
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QPushButton, 
    QLineEdit, QLabel, QGridLayout, QScrollArea, QCheckBox, QFileDialog,
    QFrame
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from pathlib import Path
import pickle


class EnhancedMeshViewer(QWidget):
    # Configuration variables for default directories

    DEFAULT_MESH_DIR = os.path.expanduser("~/dev/riboxyz/wenjun_data_meshes")
    DEFAULT_LANDMARKS_DIR = os.path.expanduser("~/dev/riboxyz/wenjun_data_landmarks")
    DEFAULT_RIBOSOME_PROFILES_DIR = os.path.expanduser("~/dev/riboxyz/wenjun_data_profiles")
    
    # Landmark visualization parameters
    LANDMARK_RADIUS = 4  # Reduced by half from the original size (assuming original was ~5)
    
    def __init__(self):
        super().__init__()
        self.meshes = {}
        self.landmarks = {}
        self.ribosome_profiles = {}
        self.plotter = None
        self.grid_plotters = []
        self.current_grid_page = 0
        self.meshes_per_page = 20  # 3x3 grid
        self.show_landmarks = True
        
        # Set default directories from config variables
        self.mesh_dir = self.DEFAULT_MESH_DIR
        self.landmarks_dir = self.DEFAULT_LANDMARKS_DIR
        self.ribosome_profiles_dir = self.DEFAULT_RIBOSOME_PROFILES_DIR
        
        # Create a monospace font for all annotations
        self.mono_font = QFont("Courier New")
        self.mono_font.setStyleHint(QFont.Monospace)
        
        self.init_ui()
        
    def init_ui(self):
        self.setWindowTitle('Enhanced PLY Mesh Viewer')
        self.setGeometry(100, 100, 1200, 800)
        
        main_layout = QVBoxLayout()
        
        # Control panel - top row
        control_panel_top = QHBoxLayout()
        
        # Mesh directory selection
        self.mesh_dir_label = QLabel('Mesh Directory:')
        self.mesh_dir_label.setFont(self.mono_font)
        self.mesh_dir_input = QLineEdit(self.mesh_dir)
        self.mesh_dir_browse = QPushButton('Browse...')
        self.mesh_dir_browse.clicked.connect(self.browse_mesh_dir)
        
        # Landmarks directory selection
        self.landmarks_dir_label = QLabel('Landmarks Directory:')
        self.landmarks_dir_label.setFont(self.mono_font)
        self.landmarks_dir_input = QLineEdit(self.landmarks_dir)
        self.landmarks_dir_browse = QPushButton('Browse...')
        self.landmarks_dir_browse.clicked.connect(self.browse_landmarks_dir)
        
        # Ribosome profiles directory selection
        self.ribosome_profiles_dir_label = QLabel('Ribosome Profiles:')
        self.ribosome_profiles_dir_label.setFont(self.mono_font)
        self.ribosome_profiles_dir_input = QLineEdit(self.ribosome_profiles_dir)
        self.ribosome_profiles_dir_browse = QPushButton('Browse...')
        self.ribosome_profiles_dir_browse.clicked.connect(self.browse_ribosome_profiles_dir)
        
        # Add widgets to top control panel
        control_panel_top.addWidget(self.mesh_dir_label)
        control_panel_top.addWidget(self.mesh_dir_input)
        control_panel_top.addWidget(self.mesh_dir_browse)
        control_panel_top.addWidget(self.landmarks_dir_label)
        control_panel_top.addWidget(self.landmarks_dir_input)
        control_panel_top.addWidget(self.landmarks_dir_browse)
        control_panel_top.addWidget(self.ribosome_profiles_dir_label)
        control_panel_top.addWidget(self.ribosome_profiles_dir_input)
        control_panel_top.addWidget(self.ribosome_profiles_dir_browse)
        
        # Control panel - bottom row
        control_panel_bottom = QHBoxLayout()
        
        # Load button
        self.load_btn = QPushButton('Load Data')
        self.load_btn.clicked.connect(self.load_data)
        
        # Input for specific mesh ID
        self.id_label = QLabel('Mesh ID:')
        self.id_label.setFont(self.mono_font)
        self.id_input = QLineEdit()
        self.view_btn = QPushButton('View Single Mesh')
        self.view_btn.clicked.connect(self.view_single_mesh)
        
        # Grid navigation
        self.prev_btn = QPushButton('Previous Page')
        self.prev_btn.clicked.connect(self.prev_page)
        self.next_btn = QPushButton('Next Page')
        self.next_btn.clicked.connect(self.next_page)
        self.page_label = QLabel('Page: 1')
        self.page_label.setFont(self.mono_font)
        
        # Toggle landmarks
        self.landmarks_check = QCheckBox('Show Landmarks')
        self.landmarks_check.setChecked(True)
        self.landmarks_check.stateChanged.connect(self.toggle_landmarks)
        
        # Add widgets to bottom control panel
        control_panel_bottom.addWidget(self.load_btn)
        control_panel_bottom.addWidget(self.id_label)
        control_panel_bottom.addWidget(self.id_input)
        control_panel_bottom.addWidget(self.view_btn)
        control_panel_bottom.addWidget(self.prev_btn)
        control_panel_bottom.addWidget(self.page_label)
        control_panel_bottom.addWidget(self.next_btn)
        control_panel_bottom.addWidget(self.landmarks_check)
        
        main_layout.addLayout(control_panel_top)
        main_layout.addLayout(control_panel_bottom)
        
        # Grid layout for mesh thumbnails
        self.grid_container = QWidget()
        self.grid_layout = QGridLayout(self.grid_container)
        
        # Scroll area for grid
        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_area.setWidget(self.grid_container)
        
        main_layout.addWidget(scroll_area)
        
        self.setLayout(main_layout)
    
    def browse_mesh_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Mesh Directory")
        if dir_path:
            self.mesh_dir_input.setText(dir_path)
            self.mesh_dir = dir_path
    
    def browse_landmarks_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Landmarks Directory")
        if dir_path:
            self.landmarks_dir_input.setText(dir_path)
            self.landmarks_dir = dir_path
    
    def browse_ribosome_profiles_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Ribosome Profiles Directory")
        if dir_path:
            self.ribosome_profiles_dir_input.setText(dir_path)
            self.ribosome_profiles_dir = dir_path
    
    def toggle_landmarks(self, state):
        self.show_landmarks = (state == Qt.Checked)
        self.update_grid()
    
    def load_data(self):
        # Get directory paths from input fields
        self.mesh_dir = self.mesh_dir_input.text()
        self.landmarks_dir = self.landmarks_dir_input.text()
        self.ribosome_profiles_dir = self.ribosome_profiles_dir_input.text()
        
        if not os.path.exists(self.mesh_dir):
            print(f"Mesh directory not found: {self.mesh_dir}")
            return
        
        if not os.path.exists(self.landmarks_dir):
            print(f"Landmarks directory not found: {self.landmarks_dir}")
            return
        
        if not os.path.exists(self.ribosome_profiles_dir):
            print(f"Ribosome profiles directory not found: {self.ribosome_profiles_dir}")
            return
        
        # Clear existing data
        self.meshes = {}
        self.landmarks = {}
        self.ribosome_profiles = {}
        for plotter in self.grid_plotters:
            plotter.close()
        self.grid_plotters = []
        
        # Clear grid layout
        for i in reversed(range(self.grid_layout.count())):
            widget = self.grid_layout.itemAt(i).widget()
            if widget is not None:
                widget.setParent(None)
        
        # Load meshes
        mesh_files = glob.glob(os.path.join(self.mesh_dir, '*_NPET_MESH.ply'))
        for file_path in mesh_files:
            mesh_id = os.path.basename(file_path).split('_')[0]
            self.meshes[mesh_id] = file_path
        
        # Load landmarks
        landmark_files = glob.glob(os.path.join(self.landmarks_dir, '*_landmarks_pts.json'))
        for file_path in landmark_files:
            mesh_id = os.path.basename(file_path).split('_')[0]
            if mesh_id in self.meshes:  # Only load landmarks for meshes we have
                try:
                    with open(file_path, 'r') as f:
                        landmark_data = json.load(f)
                    self.landmarks[mesh_id] = landmark_data
                except Exception as e:
                    print(f"Error loading landmarks for {mesh_id}: {e}")
        
        # Load ribosome profiles
        profile_files = glob.glob(os.path.join(self.ribosome_profiles_dir, '*.pkl'))
        for file_path in profile_files:
            try:
                mesh_id = os.path.basename(file_path).split('.')[0]
                if mesh_id in self.meshes:  # Only load profiles for meshes we have
                    with open(file_path, 'rb') as f:
                        profile_data = pickle.load(f)
                    self.ribosome_profiles[mesh_id] = profile_data
            except Exception as e:
                print(f"Error loading ribosome profile for {mesh_id}: {e}")
        
        print(f"Loaded {len(self.meshes)} meshes, {len(self.landmarks)} landmark files, and {len(self.ribosome_profiles)} ribosome profiles.")
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
            # Changed show_edges to False to remove wireframe
            plotter.add_mesh(mesh, color='lightblue', show_edges=False, opacity=0.7)
            
            # Add landmarks if enabled and available
            if self.show_landmarks and mesh_id in self.landmarks:
                try:
                    self.add_landmarks_to_plotter(plotter, mesh_id)
                except Exception as e:
                    print(f"Error adding landmarks for {mesh_id}: {e}")
                    
            plotter.add_text(mesh_id, position='upper_left', font_size=10, font='courier')
            plotter.reset_camera()
            
            # Add plotter to grid
            self.grid_layout.addWidget(plotter, row, col)
            self.grid_plotters.append(plotter)
            
            # Add a frame to contain the labels
            id_frame = QFrame()
            id_frame_layout = QVBoxLayout(id_frame)
            id_frame_layout.setContentsMargins(0, 0, 0, 0)
            id_frame_layout.setSpacing(0)
            
            # Add ID label
            id_label = QLabel(f"ID: {mesh_id}")
            id_label.setFont(self.mono_font)
            id_label.setAlignment(Qt.AlignCenter)
            id_label.setStyleSheet("background-color: rgba(200, 200, 200, 150); padding: 2px;")
            id_frame_layout.addWidget(id_label)
            
            # Add taxonomic and mitochondrial info if available
            if mesh_id in self.ribosome_profiles:
                profile = self.ribosome_profiles[mesh_id]
                
                # Extract taxonomic name (first organism from the list)
                taxonomy = "Unknown"
                if hasattr(profile, 'src_organism_names') and profile.src_organism_names:
                    full_name = profile.src_organism_names[0]
                    # Get genus and species for cleaner display
                    name_parts = full_name.split()
                    if len(name_parts) >= 2:
                        taxonomy = f"{name_parts[0][0]}. {name_parts[1]}"
                    else:
                        taxonomy = full_name[:20]
                
                # Extract mitochondrial status
                is_mito = False
                if hasattr(profile, 'mitochondrial'):
                    is_mito = profile.mitochondrial
                    
                # Color code based on mitochondrial status
                bg_color = "rgba(144, 238, 144, 180)" if is_mito else "rgba(173, 216, 230, 180)"
                mito_label = "MITO" if is_mito else "CYTO"
                
                # Create info label with color coding
                info_label = QLabel(f"{taxonomy} | {mito_label}")
                info_label.setFont(self.mono_font)
                info_label.setAlignment(Qt.AlignCenter)
                info_label.setStyleSheet(f"background-color: {bg_color}; padding: 2px; border-radius: 2px;")
                id_frame_layout.addWidget(info_label)
                
                # Add a second label with resolution if available
                if hasattr(profile, 'resolution'):
                    res_label = QLabel(f"Res: {profile.resolution:.1f}Å")
                    res_label.setFont(self.mono_font)
                    res_label.setAlignment(Qt.AlignCenter)
                    res_label.setStyleSheet("background-color: rgba(200, 200, 200, 150); padding: 2px;")
                    id_frame_layout.addWidget(res_label)
            
            # Add the frame to the grid layout in the cell below the plotter
            self.grid_layout.addWidget(id_frame, row + 1, col)
        
        # Update page label
        total_pages = max(1, int(np.ceil(len(self.meshes) / self.meshes_per_page)))
        self.page_label.setText(f'Page: {self.current_grid_page + 1}/{total_pages}')
        
        # Enable/disable navigation buttons
        self.prev_btn.setEnabled(self.current_grid_page > 0)
        self.next_btn.setEnabled(self.current_grid_page < total_pages - 1)
    
    def add_landmarks_to_plotter(self, plotter, mesh_id):
        """Add pre-generated landmarks to the plotter with reduced size"""
        landmarks_data = self.landmarks[mesh_id]
        
        # Add PTC as a yellow sphere (half the size)
        if landmarks_data.get("ptc"):
            ptc_location = np.array(landmarks_data["ptc"])
            plotter.add_mesh(pv.Sphere(center=ptc_location, radius=self.LANDMARK_RADIUS), color='yellow')
        
        # Add constriction site as a red sphere (half the size)
        if landmarks_data.get("constriction"):
            constriction_location = np.array(landmarks_data["constriction"])
            plotter.add_mesh(pv.Sphere(center=constriction_location, radius=self.LANDMARK_RADIUS), color='red')
        
        # Add uL4 residues as blue points (half the size)
        if landmarks_data.get("uL4"):
            ul4_points = np.array(landmarks_data["uL4"])
            if len(ul4_points) > 0:
                plotter.add_mesh(pv.PolyData(ul4_points), color='blue', 
                                point_size=5,  # Reduced from 10
                                render_points_as_spheres=True)
        
        # Add uL22 residues as orange points (half the size)
        if landmarks_data.get("uL22"):
            ul22_points = np.array(landmarks_data["uL22"])
            if len(ul22_points) > 0:
                plotter.add_mesh(pv.PolyData(ul22_points), color='orange', 
                                point_size=5,  # Reduced from 10
                                render_points_as_spheres=True)
    
    def prev_page(self):
        if self.current_grid_page > 0:
            self.current_grid_page -= 1
            self.update_grid()
    
    def next_page(self):
        total_pages = int(np.ceil(len(self.meshes) / self.meshes_per_page))
        if self.current_grid_page < total_pages - 1:
            self.current_grid_page += 1
            self.update_grid()
    
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
        
        # Create new plotter window that won't kill the application when closed
        self.plotter = pv.Plotter(off_screen=False)
        
        # Load and display mesh
        mesh = pv.read(self.meshes[mesh_id])
        # Changed show_edges to False to remove wireframe
        self.plotter.add_mesh(mesh, color='lightblue', show_edges=False, opacity=0.7)
        
        # Add landmarks if enabled and available
        if self.show_landmarks and mesh_id in self.landmarks:
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
        
        # Add ribosome profile information if available
        if mesh_id in self.ribosome_profiles:
            profile = self.ribosome_profiles[mesh_id]
            
            # Add taxonomic information
            if hasattr(profile, 'src_organism_names') and profile.src_organism_names:
                org_name = profile.src_organism_names[0]
                info_text.append(f"Organism: {org_name}")
                
                # Add taxonomy ID if available
                if hasattr(profile, 'src_organism_ids') and profile.src_organism_ids:
                    taxid = profile.src_organism_ids[0]
                    info_text.append(f"Taxonomy ID: {taxid}")
            
            # Add mitochondrial status with clear indicator
            if hasattr(profile, 'mitochondrial'):
                is_mito = profile.mitochondrial
                status = "MITOCHONDRIAL" if is_mito else "CYTOSOLIC"
                info_text.append(f"Type: {status}")
            
            # Add resolution with Angstrom symbol
            if hasattr(profile, 'resolution'):
                info_text.append(f"Resolution: {profile.resolution:.2f}Å")
                
            # Add experiment method if available
            if hasattr(profile, 'expMethod'):
                info_text.append(f"Method: {profile.expMethod}")
                
            # Add publication info if available
            if hasattr(profile, 'citation_title') and profile.citation_title:
                info_text.append(f"Publication: {profile.citation_title}")
            
            # Add deposition date if available
            if hasattr(profile, 'deposition_date') and profile.deposition_date:
                info_text.append(f"Deposited: {profile.deposition_date}")
        
        # Add landmark status information
        if mesh_id in self.landmarks:
            lm_data = self.landmarks[mesh_id]
            lm_info = []
            if lm_data.get("ptc"):
                lm_info.append("PTC: ✓")
            else:
                lm_info.append("PTC: ✗")
                
            if lm_data.get("constriction"):
                lm_info.append("Constriction: ✓")
            else:
                lm_info.append("Constriction: ✗")
                
            if lm_data.get("uL4") and len(lm_data["uL4"]) > 0:
                lm_info.append(f"uL4: ✓ ({len(lm_data['uL4'])} points)")
            else:
                lm_info.append("uL4: ✗")
                
            if lm_data.get("uL22") and len(lm_data["uL22"]) > 0:
                lm_info.append(f"uL22: ✓ ({len(lm_data['uL22'])} points)")
            else:
                lm_info.append("uL22: ✗")
                
            info_text.extend(lm_info)
        else:
            info_text.append("No landmark data available")
            
        self.plotter.add_text('\n'.join(info_text), position='upper_left', font_size=12, font='courier')
        
        # Show the plotter in a non-blocking way that won't kill the app when closed
        self.plotter.show(title=f"Mesh Viewer - {mesh_id}", auto_close=False)


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
            print("pip install pyvista pyvistaqt numpy PyQt5")