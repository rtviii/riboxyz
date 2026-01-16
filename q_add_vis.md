Hey, so i have this pipeline orchestrator for processing a ribosome exit tunnel and creating a mesh from it, but it involves a few stages and i want to be able to register different visualizations for each stage + some visualization that don't depend on any stage, but can just be toggled on and off. Let me describe to you the layout of my project.

Here is the project directory (most of the stuff that matters is under the "npet" folder):
```
ᢹ saeta.rtviii[ dev/riboxyz ]  tree -I 'node_modules|venv|*__pycache__*|profiles|cache|debug_output|*json|*npy|*ply|*fasta|*csv|assets_*|staticfiles|api|assets' -L 5                                                             [meshing_methods_refactor]
.
├── __scripts
│   ├── gettrna.py
│   ├── ligclass_distill_composite.py
│   ├── ligclass_distill.py
│   ├── ligclass.py
│   ├── maps_acquire.py
│   ├── masif_plugin.py
│   ├── merge_lig_info.py
│   ├── mesh_grid.py
│   ├── move_classification.sh
│   ├── move_tunnels.sh
│   ├── narstructs_tally
│   ├── npy_to_pdb.py
│   ├── process_prediction_7k00.py
│   ├── pymol_visualtion.py
│   ├── rcsb_cumulative_entries_barplot.py
│   ├── reclassify.py
│   ├── ribxz_chimerax
│   │   ├── build
│   │   │   ├── bdist.macosx-10.9-universal2
│   │   │   └── lib
│   │   │       └── ribxz_chimerax
│   │   ├── bundle_info.xml
│   │   ├── dist
│   │   │   └── ribxz_chimerax-0.1-py3-none-any.whl
│   │   ├── ribxz_chimerax.egg-info
│   │   │   ├── dependency_links.txt
│   │   │   ├── PKG-INFO
│   │   │   ├── requires.txt
│   │   │   ├── SOURCES.txt
│   │   │   └── top_level.txt
│   │   └── src
│   │       ├── __init__.py
│   │       ├── cmd_loci.py
│   │       ├── cmd_polymers.py
│   │       ├── cmd_registry.py
│   │       └── io.py
│   ├── uniprot_seeds_query_record.py
│   └── williamson_assembly.py
├── docker-compose.yml
├── Dockerfile-django
├── docs.md
├── kingdom_analysis.py
├── neo4j_ribosome
│   ├── __archive
│   │   ├── cypher_ops
│   │   │   ├── cypher_exec
│   │   │   ├── neo4j_commit_structure.sh
│   │   │   └── neo4j_seed_db_ontology.sh
│   │   └── riboxyz_seed_data
│   ├── __init__.py
│   ├── cypher
│   │   └── list_filtered_structs.cypher
│   ├── db_driver.py
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── node_ligand.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_protein.py
│   ├── node_rna.py
│   └── node_structure.py
├── neo4j.conf
├── neo4j.old.conf
├── notes
│   ├── binding_affinities.md
│   ├── docs
│   │   ├── architecture_environment_configs.md
│   │   └── update.md
│   ├── factors.md
│   ├── functional_sites.md
│   ├── general_directions.md
│   ├── hmm-based-classification.md
│   └── pymol.md
├── npet_bulk_viewer.py
├── npet_inspect.py
├── npet_orchestrator.py
├── npet_plan.md
├── pipeline_manager.py
├── q_add_vis.md
├── q_orchestrator.md
├── q.md
├── ribctl
│   ├── __init__.py
│   ├── asset_manager
│   │   ├── asset_manager.py
│   │   ├── asset_registry.py
│   │   ├── asset_types.py
│   │   ├── doc.md
│   │   └── parallel_acquisition.py
│   ├── etl
│   │   ├── __init__.py
│   │   ├── etl_collector.py
│   │   └── gql_querystrings.py
│   ├── global_ops.py
│   ├── lib
│   │   ├── __libseq.py
│   │   ├── chimerax
│   │   │   ├── _cmd_ribrepr.py
│   │   │   ├── cmd_chainsplitter.py
│   │   │   ├── cmd_ligvis.py
│   │   │   ├── cmd_ribetl.py
│   │   │   ├── cmd_ribmovie.py
│   │   │   ├── cmd_ribrepr.py
│   │   │   ├── cmds_all.py
│   │   │   ├── ffmpeg_convert.sh
│   │   │   ├── ffmpeg_firstframe.sh
│   │   │   ├── gen_movies.py
│   │   │   ├── loop_ligvis.py
│   │   │   ├── loop_movies.py
│   │   │   ├── loop_split_chains.py
│   │   │   ├── notes.md
│   │   │   ├── produce_gif.py
│   │   │   └── thumbnails_from_ribetl.sh
│   │   ├── enumunion.py
│   │   ├── info.py
│   │   ├── landmarks
│   │   │   ├── __ptc_via_doris.py
│   │   │   ├── constriction_site.py
│   │   │   ├── notes.md
│   │   │   ├── ptc_via_trna.py
│   │   │   └── rrna_helices
│   │   │       ├── convert.py
│   │   │       └── rrna_helices.py
│   │   ├── libbsite.py
│   │   ├── libhmm.py
│   │   ├── libmsa.py
│   │   ├── libseq.py
│   │   ├── libtax.py
│   │   ├── npet
│   │   │   ├── _dbscan_pairs
│   │   │   ├── alphalib.py
│   │   │   ├── kdtree_approach.py
│   │   │   ├── landmark_ptcloud_collector.py
│   │   │   ├── npet_pipeline.py
│   │   │   ├── pipeline
│   │   │   │   ├── base_stage.py
│   │   │   │   ├── clustering_stage.py
│   │   │   │   ├── constriction_identification_stage.py
│   │   │   │   ├── entity_filtering_stage.py
│   │   │   │   ├── exterior_mesh_stage.py
│   │   │   │   ├── landmark_identification_stage.py
│   │   │   │   ├── logs
│   │   │   │   ├── mesh_reconstruction_stage.py
│   │   │   │   ├── normal_estimation_stage.py
│   │   │   │   ├── point_cloud_processing_stage.py
│   │   │   │   ├── ptc_identification_stage.py
│   │   │   │   ├── refinement_stage.py
│   │   │   │   ├── setup_stage.py
│   │   │   │   ├── surface_extraction_stage.py
│   │   │   │   └── validation_stage.py
│   │   │   ├── pipeline_status_tracker.py
│   │   │   ├── tunnel_asset_manager.py
│   │   │   └── various_visualization.py
│   │   ├── nsearch_gemmi.py
│   │   ├── ribosome_types
│   │   ├── schema
│   │   │   ├── __init__.py
│   │   │   ├── primitives.py
│   │   │   ├── types_binding_site.py
│   │   │   └── types_ribosome.py
│   │   ├── seq_project_many_to_one.py
│   │   ├── thumbnail.py
│   │   ├── types
│   │   │   └── polymer
│   │   │       ├── __init__.py
│   │   │       ├── base.py
│   │   │       ├── hierarchies.py
│   │   │       └── types.py
│   │   └── utils.py
│   ├── logger_config.py
│   ├── logs
│   │   ├── etl.log
│   │   └── loggers.py
│   ├── ribd.py
│   └── ribosome_ops.py
├── ribxz_logo_black.png
├── taxdump.tar.gz
└── visualization_library.py

32 directories, 148 files
```
The thing im building constitutes this "npet_orchestraotr" and "pipeline_manager" files more or less, but they call stuff in the "npet" library, which is in the "ribctl" folder. The pipeline is broken down into stages and each stage relies on some piece of computation takes in some artifactsa and produces some artifacts, which, along with the primary data (the ribosome structures and such) are stored in the general RIBETL_DATA dir. The artifacts path structure is typed out in the "npet/tunnel_asset_manager.py" file. Each stage is defined as a class in the "npet/pipeline" folder, and they all inherit from a "BaseStage" class that defines the general interface for a stage. Each stage has a "run" method that does the computation and produces the artifacts.:
```
import os

from ribctl import ASSETS_PATH, RIBETL_DATA, RIBXZ_TEMP_FILES
from ribctl.ribosome_ops import RibosomeOps


class TunnelMeshAssetsManager:

    rcsb_id: str

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id.upper()
        RO = RibosomeOps(rcsb_id)
        self.structpath = RO.assets.paths.cif
        pass

    @property
    def cif_struct(self):
        return self.structpath

    @property
    def tunnel_pcd_normal_estimated(self):
        return os.path.join( RIBXZ_TEMP_FILES, "{}_tunnel_pcd_normal_estimated.ply".format(self.rcsb_id) )

    @property
    def tunnel_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_tunnel_half_poisson_recon.ply".format(self.rcsb_id)
        )

    @property
    def ashape_half_mesh(self):
        return os.path.join(
            self.structpath, "{}_half_ashape_watertight.ply".format(self.rcsb_id)
        )

    # These are just to visualize lamps trajectories/pointclouds as pdb files.

    @property
    def lammps_traj_tunnel(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj".format(self.rcsb_id)
        )

    @property
    def lammps_traj_tunnel_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_tunnel.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def lammps_traj_ashape_as_pdb(self):
        return os.path.join(
            self.structpath, "{}_ashape.lammpstraj.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_PDB_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.pdb".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_mmcif_largest(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.largest.cif".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_noise(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.noise.xyz".format(self.rcsb_id)
        )

    @property
    def dbscan_clusters_xyz_refined(self):
        return os.path.join(
            self.structpath, "{}_dbscan_clusters.refined.xyz".format(self.rcsb_id)
        )
```


Basically what i want to do right now if you can is take one stage (PTC identification) and implement a custom visualizaiton for it, but as a proof of concept for a general visualziation regsitration system. Simultaneously we want to add a button the ui that lets us toggle on and off the general riboosme structure visualization -- that is the ui, no matter which stage we are on, should be able to grab the ribosome cif, parse its resiudues cneters and add those as balls to the viewer.  Can we try this?

npet_orchestrator.py:
```
#!/usr/bin/env python
"""
NPET Pipeline Orchestrator
Interactive GUI for running and debugging the NPET mesh generation pipeline.
"""

import sys
sys.path.append("/home/rtviii/dev/riboxyz")
import os
import json
import traceback
from pathlib import Path
from typing import Dict, Any, Optional, List, Callable
from enum import Enum
import pickle
from datetime import datetime
from dataclasses import dataclass, field # <-- Added

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, 
    QPushButton, QLabel, QTextEdit, QSplitter, QScrollArea, QFrame,
    QLineEdit, QGroupBox, QGridLayout, QSpinBox, QDoubleSpinBox,
    QComboBox, QFileDialog, QProgressBar, QMessageBox, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont

# --- IMPORT THE NEW DEFINITIONS ---
# These must now exist in your pipeline_manager.py file
from pipeline_manager import PipelineManager, StageDefinition, VisualizationSpec

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
    (Refactored from visualize_mesh)
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
    background_mesh: Optional[pv.PolyData] = None, # <-- ADD THIS ARGUMENT
    point_cloud_data: Optional[np.ndarray] = None,
    ptc_coord: Optional[np.ndarray] = None,
    constriction_coord: Optional[np.ndarray] = None,
):
    """
    Visualizes landmarks and/or a point cloud,
    optionally on top of a background mesh.
    """
    
    # --- ADD THIS BLOCK ---
    # 1. Add the background mesh if provided
    if background_mesh is not None:
        plotter.add_mesh(
            background_mesh,
            style='surface',
            color='white',
            opacity=0.1, # <-- Make it faint
            show_edges=False
        )
    # --- END OF BLOCK ---

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
# -------------------------------------------------------------------
## --- CORE GUI CLASSES ---
# -------------------------------------------------------------------

class StageStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"

class StageExecutor(QThread):
    """
    Worker thread for running a pipeline stage.
    """
    finished = pyqtSignal(str, bool, object, str) 
    
    def __init__(self, manager: PipelineManager, stage_name: str, params: Dict[str, Any]):
        super().__init__()
        self.manager = manager
        self.stage_name = stage_name
        self.params = params
        
    def run(self):
        try:
            result = self.manager.run_stage(self.stage_name, self.params)
            self.finished.emit(self.stage_name, True, result, "Stage completed successfully")
        except Exception as e:
            error_msg = f"Stage failed: {str(e)}\n{traceback.format_exc()}"
            self.finished.emit(self.stage_name, False, None, error_msg)

# --- StageWidget (Unchanged from your file) ---
class StageWidget(QFrame):
    """Widget representing a single pipeline stage"""
    
    run_requested = pyqtSignal(StageDefinition)
    view_artifacts = pyqtSignal(StageDefinition)
    
    def __init__(self, config: StageDefinition, parent=None):
        super().__init__(parent)
        self.config = config
        self.status = StageStatus.PENDING
        self.params = config.default_params.copy()
        
        self.artifacts: Dict[str, Path] = {} 
        self.execution_time = None
        self.setup_ui()
        self.update_status_display()
        
    def setup_ui(self):
        self.setFrameStyle(QFrame.Box | QFrame.Raised)
        self.setLineWidth(2)
        
        layout = QVBoxLayout(self)
        
        # Header
        header_layout = QHBoxLayout()
        self.status_label = QLabel("●")
        self.status_label.setFont(QFont("Arial", 16))
        header_layout.addWidget(self.status_label)
        name_label = QLabel(self.config.name)
        name_label.setFont(QFont("Arial", 12, QFont.Bold))
        header_layout.addWidget(name_label)
        header_layout.addStretch()
        self.time_label = QLabel("")
        self.time_label.setFont(QFont("Courier New", 9))
        header_layout.addWidget(self.time_label)
        layout.addLayout(header_layout)
        
        # Description
        desc_label = QLabel(self.config.description)
        desc_label.setWordWrap(True)
        desc_label.setStyleSheet("color: gray; font-size: 10px;")
        layout.addWidget(desc_label)
        
        # Parameters
        params_group = QGroupBox("Parameters")
        params_layout = QGridLayout()
        self.param_widgets = {}
        row = 0
        for param_name, param_value in self.params.items():
            label = QLabel(f"{param_name}:")
            params_layout.addWidget(label, row, 0)
            
            if isinstance(param_value, bool):
                widget = QComboBox()
                widget.addItems(["True", "False"])
                widget.setCurrentText(str(param_value))
            elif isinstance(param_value, int):
                widget = QSpinBox()
                widget.setRange(-10000, 10000)
                widget.setValue(param_value)
            elif isinstance(param_value, float):
                widget = QDoubleSpinBox()
                widget.setRange(-1000.0, 1000.0)
                widget.setSingleStep(0.1)
                widget.setDecimals(3)
                widget.setValue(param_value)
            else:
                widget = QLineEdit(str(param_value))
            
            self.param_widgets[param_name] = widget
            params_layout.addWidget(widget, row, 1)
            row += 1
        
        params_group.setLayout(params_layout)
        layout.addWidget(params_group)
        
        # Artifacts info
        self.artifacts_label = QLabel("Artifacts: 0")
        self.artifacts_label.setFont(QFont("Courier New", 9))
        layout.addWidget(self.artifacts_label)
        
        # Action buttons
        button_layout = QHBoxLayout()
        self.run_btn = QPushButton("Run Stage")
        self.run_btn.clicked.connect(lambda: self.run_requested.emit(self.config))
        button_layout.addWidget(self.run_btn)
        
        self.view_btn = QPushButton("View Artifacts")
        self.view_btn.setEnabled(False)
        self.view_btn.clicked.connect(lambda: self.view_artifacts.emit(self.config))
        button_layout.addWidget(self.view_btn)
        
        layout.addLayout(button_layout)

    def update_status_display(self):
        colors = {
            StageStatus.PENDING: "gray",
            StageStatus.RUNNING: "orange",
            StageStatus.SUCCESS: "green",
            StageStatus.FAILED: "red",
            StageStatus.SKIPPED: "lightgray"
        }
        color = colors.get(self.status, "gray")
        self.status_label.setStyleSheet(f"color: {color};")
        
        self.run_btn.setEnabled(self.status != StageStatus.RUNNING)
        self.view_btn.setEnabled(bool(self.artifacts) and self.status != StageStatus.RUNNING)
        
        self.artifacts_label.setText(f"Artifacts: {len(self.artifacts)}")
        if self.execution_time:
            self.time_label.setText(f"{self.execution_time:.2f}s")
        else:
            self.time_label.setText("")

    def set_status(self, status: StageStatus, time: Optional[float] = None):
        self.status = status
        if time is not None:
            self.execution_time = time
        self.update_status_display()
        
    def get_params(self) -> Dict[str, Any]:
        params = {}
        for name, widget in self.param_widgets.items():
            if isinstance(widget, QComboBox):
                params[name] = widget.currentText() == "True"
            elif isinstance(widget, QSpinBox):
                params[name] = widget.value()
            elif isinstance(widget, QDoubleSpinBox):
                params[name] = widget.value()
            else: # QLineEdit
                val = widget.text()
                try:
                    params[name] = int(val)
                except ValueError:
                    try:
                        params[name] = float(val)
                    except ValueError:
                        params[name] = val
        return params

    def set_artifacts(self, artifacts: Dict[str, Path]):
        self.artifacts = artifacts
        self.update_status_display()

# --- ArtifactViewer (HEAVILY REFACTORED) ---
class ArtifactViewer(QWidget):
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.manager: Optional[PipelineManager] = None
        self.stage_config: Optional[StageDefinition] = None
        self.current_params: Dict[str, Any] = {}
        # This will hold ALL discoverable artifacts, not just for the current stage
        self.all_artifacts: Dict[str, Path] = {} 
        
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        info_group = QGroupBox("Visualization Info")
        info_layout = QVBoxLayout()
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMinimumHeight(150)
        self.info_text.setFont(QFont("Courier New", 9))
        info_layout.addWidget(self.info_text)
        info_group.setLayout(info_layout)
        layout.addWidget(info_group)
        
        selector_layout = QGridLayout()
        
        selector_layout.addWidget(QLabel("Visualization:"), 0, 0)
        self.viz_select_combo = QComboBox()
        self.viz_select_combo.currentIndexChanged.connect(self.on_viz_selected)
        selector_layout.addWidget(self.viz_select_combo, 0, 1)
        
        self.raw_file_label = QLabel("Select File:")
        self.raw_file_label.hide()
        selector_layout.addWidget(self.raw_file_label, 1, 0)
        
        self.raw_file_combo = QComboBox()
        self.raw_file_combo.currentIndexChanged.connect(self.load_raw_artifact)
        self.raw_file_combo.hide()
        selector_layout.addWidget(self.raw_file_combo, 1, 1)

        self.export_btn = QPushButton("Export...")
        self.export_btn.clicked.connect(self.export_artifact)
        self.export_btn.hide()
        selector_layout.addWidget(self.export_btn, 1, 2)
        
        layout.addLayout(selector_layout)
        
        self.plotter_frame = QFrame()
        self.plotter_frame.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        plotter_layout = QVBoxLayout(self.plotter_frame)
        
        self.plotter = QtInteractor(self.plotter_frame)
        self.plotter.enable_trackball_style()
        plotter_layout.addWidget(self.plotter)  
        
        layout.addWidget(self.plotter_frame, 1)
        
    def set_stage_context(self, 
                          stage_config: StageDefinition, 
                          manager: PipelineManager, 
                          current_params: Dict[str, Any]):
        """
        Called by the main window to set the *context* for the viewer.
        This populates the dropdowns.
        """
        self.stage_config = stage_config
        self.manager = manager
        self.current_params = current_params
        
        # --- *** ONE-LINE MODIFICATION *** ---
        # This is CRITICAL. It ensures that self.manager.context is
        # populated from disk artifacts *before* we try to use it.
        try:
            self.manager._load_context_for_stage(stage_config)
        except Exception as e:
            print(f"Warning: Failed to pre-load context: {e}")
            # Don't crash, just log it. The viz function will fail later.
        
        # Build a map of ALL artifacts from ALL stages
        self.all_artifacts.clear()
        for stage in self.manager.stages:
            self.all_artifacts.update(self.manager.get_stage_artifacts(stage.name))

        # Get artifacts for *this* stage (for the raw inspector)
        current_stage_artifacts = self.manager.get_stage_artifacts(stage_config.name)
        
        self.viz_select_combo.blockSignals(True)
        self.raw_file_combo.blockSignals(True)
        
        self.viz_select_combo.clear()
        self.raw_file_combo.clear()
        
        self.info_text.setPlainText(f"Stage: {stage_config.name}\n"
                                    f"Found {len(current_stage_artifacts)} artifact files for this stage.")

        # 1. Add registered visualizations
        for viz_spec in stage_config.visualizations:
            self.viz_select_combo.addItem(viz_spec.name, userData=viz_spec)
            
        # 2. Add the "Raw File Inspector" as a fallback
        self.viz_select_combo.addItem("Inspect Raw Artifact File", userData="RAW_FILE_INSPECTOR")

        # 3. Populate the (hidden) raw file dropdown
        for name, path in current_stage_artifacts.items():
            if path.exists():
                self.raw_file_combo.addItem(f"{name} ({path.name})", userData=path)
        
        self.viz_select_combo.blockSignals(False)
        self.raw_file_combo.blockSignals(False)
        
        self.on_viz_selected(0)

    def on_viz_selected(self, index: int):
        """Handles switching between different visualizations."""
        if index < 0:
            return
            
        viz_data = self.viz_select_combo.itemData(index)
        
        if viz_data == "RAW_FILE_INSPECTOR":
            self.raw_file_label.show()
            self.raw_file_combo.show()
            self.export_btn.show()
            self.load_raw_artifact(self.raw_file_combo.currentIndex())
            
        elif isinstance(viz_data, VisualizationSpec):
            self.raw_file_label.hide()
            self.raw_file_combo.hide()
            self.export_btn.hide()
            self.execute_visualization(viz_data)

    def load_data_from_path(self, path: Path) -> Any:
        """Helper to load common artifact types from disk."""
        ext = path.suffix.lower()
        try:
            if ext == '.npy':
                return np.load(path, allow_pickle=True)
            elif ext in ['.ply', '.stl', '.obj', '.vtk']:
                return pv.read(path)
            elif ext == '.json':
                with open(path, 'r') as f:
                    return json.load(f)
            elif ext == '.cif':
                return str(path)
            else:
                return str(path)
        except Exception as e:
            print(f"Error loading data from {path}: {e}")
            raise

    def execute_visualization(self, viz_spec: VisualizationSpec):
        """Gathers data and calls the registered viz function."""
        
        self.plotter.clear()
        viz_func = viz_spec.function
        kwargs = {} 
        
        info = f"Running: {viz_spec.name}\n"
        info += f"Function: {viz_func.__name__}\n\n"
        
        try:
            # 1. Gather Artifacts
            info += "Loading Artifacts:\n"
            for kwarg_name, artifact_key in viz_spec.artifact_map.items():
                # Use the 'all_artifacts' map to find the artifact
                path = self.all_artifacts.get(artifact_key) 
                if path and path.exists():
                    data = self.load_data_from_path(path)
                    kwargs[kwarg_name] = data
                    info += f" - {kwarg_name} <-- {path.name} [OK]\n"
                else:
                    kwargs[kwarg_name] = None # Pass None to the function
                    info += f" - {kwarg_name} <-- {artifact_key} [MISSING]\n"
                    print(f"Warning: Missing artifact '{artifact_key}' at: {path}")

            # 2. Gather Context
            info += "\nLoading Context:\n"
            for kwarg_name, context_key in viz_spec.context_map.items():
                value = self.manager.context.get(context_key)
                kwargs[kwarg_name] = value
                info += f" - {kwarg_name} <-- context['{context_key}'] [OK]\n"
                
            # 3. Gather Parameters
            info += "\nLoading Parameters:\n"
            for kwarg_name, param_key in viz_spec.param_map.items():
                value = self.current_params.get(param_key)
                kwargs[kwarg_name] = value
                info += f" - {kwarg_name} <-- param['{param_key}'] [OK]\n"

            # 4. Execute the function
            self.info_text.setPlainText(info + "\nExecuting visualization...")
            
            viz_func(self.plotter, **kwargs)
            
            self.plotter.reset_camera()
            self.info_text.append("Done.")
            
        except Exception as e:
            error_msg = f"Failed to execute visualization: {viz_spec.name}\n\n"
            error_msg += f"Error: {str(e)}\n"
            error_msg += f"{traceback.format_exc()}"
            self.info_text.setPlainText(error_msg)
            self.plotter.add_text(f"Error: {e}", position='upper_left', color='red')

    def load_raw_artifact(self, index: int):
        """
        This is the *old* load_selected_artifact logic,
        used for the "Raw File Inspector".
        """
        if index < 0:
            self.plotter.clear()
            self.info_text.setPlainText("No raw files available for this stage.")
            return
            
        artifact_path: Path = self.raw_file_combo.itemData(index)
        if not artifact_path or not artifact_path.exists():
            self.info_text.setPlainText(f"\nERROR: File not found: {artifact_path}")
            return
            
        self.plotter.clear()
        self.plotter.enable_trackball_style()
        
        try:
            ext = artifact_path.suffix.lower()
            
            if ext in ['.ply', '.stl', '.obj', '.vtk']:
                mesh = pv.read(artifact_path)
                self.plotter.add_mesh(mesh, color='lightblue', show_edges=True, opacity=0.8)
                info = f"Mesh: {artifact_path.name}\n"
                info += f"Points: {mesh.n_points}\n"
                info += f"Faces: {mesh.n_faces}\n"
                self.info_text.setPlainText(info)
                
            elif ext == '.npy':
                points = np.load(artifact_path)
                if len(points.shape) == 2 and points.shape[1] == 3:
                    cloud = pv.PolyData(points)
                    self.plotter.add_mesh(cloud, color='red', point_size=5,
                                          render_points_as_spheres=True)
                    info = f"Point Cloud: {artifact_path.name}\n"
                    info += f"Points: {len(points)}\n"
                    self.info_text.setPlainText(info)
                else:
                    self.info_text.setPlainText(f"Numpy array: {artifact_path.name}\nShape: {points.shape}")

            elif ext == '.cif':
                info = f"Structure File: {artifact_path.name}\n"
                info += f"Size: {artifact_path.stat().st_size / (1024*1024):.2f} MB\n\n"
                info += "Raw CIF preview not supported."
                self.info_text.setPlainText(info)
                self.plotter.add_text("CIF structure file (see info panel)", position='upper_left')

            elif ext == '.json':
                with open(artifact_path) as f:
                    data = json.load(f)
                self.info_text.setPlainText(json.dumps(data, indent=2))
                self.plotter.add_text("JSON data (see info panel)", position='upper_left')
                
            else:
                self.info_text.setPlainText(f"Cannot preview file type: {ext}")
                self.plotter.add_text(f"Cannot preview {ext}", position='upper_left')

            self.plotter.reset_camera()
            
        except Exception as e:
            self.info_text.setPlainText(f"Error loading artifact:\n{str(e)}")
            print(traceback.format_exc())

    def export_artifact(self):
        current_index = self.raw_file_combo.currentIndex()
        if current_index < 0:
            return
        
        artifact_path: Path = self.raw_file_combo.itemData(current_index)
        if not artifact_path:
            return
            
        save_path, _ = QFileDialog.getSaveFileName(
            self, "Export Artifact",
            str(Path.home() / artifact_path.name), 
            f"{artifact_path.suffix.upper()} files (*{artifact_path.suffix});;All files (*.*)"
        )
        
        if save_path:
            try:
                import shutil
                shutil.copy(artifact_path, save_path)
                self.info_text.append(f"\nExported to {save_path}")
            except Exception as e:
                self.info_text.append(f"\nExport failed: {e}")

class NPETPipelineOrchestrator(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.rcsb_id: Optional[str] = None
        self.manager: Optional[PipelineManager] = None
        self.stage_widgets: Dict[str, StageWidget] = {}
        self.current_executor: Optional[StageExecutor] = None
        
        self.setup_ui()
        self.create_menu_bar()
        
        self.setWindowTitle("NPET Pipeline Orchestrator")
        self.setGeometry(100, 100, 1400, 900)
        self.log("Orchestrator started. Enter an RCSB ID to begin.")

    def create_menu_bar(self):
        pass

    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        control_layout = QHBoxLayout()
        control_layout.addWidget(QLabel("RCSB ID:"))
        self.id_input = QLineEdit()
        self.id_input.setPlaceholderText("e.g., 7K00")
        self.id_input.setMaximumWidth(100)
        self.id_input.returnPressed.connect(self.initialize_pipeline)
        control_layout.addWidget(self.id_input)
        
        self.init_btn = QPushButton("Load / Initialize")
        self.init_btn.clicked.connect(self.initialize_pipeline)
        control_layout.addWidget(self.init_btn)
        
        self.run_all_btn = QPushButton("Run All Stages")
        self.run_all_btn.setEnabled(False)
        self.run_all_btn.clicked.connect(self.run_all_stages)
        control_layout.addWidget(self.run_all_btn)
        
        self.reset_btn = QPushButton("Reset Pipeline")
        self.reset_btn.clicked.connect(self.reset_pipeline)
        control_layout.addWidget(self.reset_btn)
        
        control_layout.addStretch()
        main_layout.addLayout(control_layout)
        
        splitter = QSplitter(Qt.Horizontal)
        
        stages_widget = QWidget()
        stages_layout_v = QVBoxLayout(stages_widget)
        stages_label = QLabel("Pipeline Stages")
        stages_label.setFont(QFont("Arial", 14, QFont.Bold))
        stages_layout_v.addWidget(stages_label)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(400)
        
        stages_container = QWidget()
        self.stages_layout = QVBoxLayout(stages_container) 
        self.stages_layout.addStretch()
        scroll.setWidget(stages_container)
        stages_layout_v.addWidget(scroll)
        splitter.addWidget(stages_widget)
        
        self.tabs = QTabWidget()
        self.artifact_viewer = ArtifactViewer()
        self.tabs.addTab(self.artifact_viewer, "Artifact Viewer")
        
        log_widget = QWidget()
        log_layout = QVBoxLayout(log_widget)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier New", 9))
        log_layout.addWidget(self.log_text)
        self.tabs.addTab(log_widget, "Execution Log")
        
        context_widget = QWidget()
        context_layout = QVBoxLayout(context_widget)
        self.context_table = QTableWidget()
        self.context_table.setColumnCount(2)
        self.context_table.setHorizontalHeaderLabels(["Key", "Value"])
        self.context_table.horizontalHeader().setStretchLastSection(True)
        context_layout.addWidget(self.context_table)
        self.tabs.addTab(context_widget, "Pipeline Context")
        
        splitter.addWidget(self.tabs)
        splitter.setSizes([450, 950])
        main_layout.addWidget(splitter)

    def log(self, message: str, level: str = "INFO"):
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted = f"[{timestamp}] {level}: {message}"
        self.log_text.append(formatted)
        print(formatted)
        
    def initialize_pipeline(self):
        rcsb_id = self.id_input.text().strip().upper()
        if not rcsb_id:
            QMessageBox.warning(self, "Error", "Please enter an RCSB ID")
            return
            
        self.log(f"Initializing pipeline for {rcsb_id}...")
        self.rcsb_id = rcsb_id
        
        try:
            self.manager = PipelineManager(self.rcsb_id)
        except Exception as e:
            self.log(f"Failed to initialize manager: {e}\n{traceback.format_exc()}", "ERROR")
            QMessageBox.critical(self, "Error", f"Failed to initialize manager:\n{e}")
            return

        while self.stages_layout.count() > 1:
            child = self.stages_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.stage_widgets.clear()

        for stage_def in self.manager.stages:
            stage_widget = StageWidget(stage_def)
            stage_widget.run_requested.connect(self.run_stage)
            stage_widget.view_artifacts.connect(self.view_stage_artifacts)
            
            self.stages_layout.insertWidget(self.stages_layout.count() - 1, stage_widget)
            self.stage_widgets[stage_def.name] = stage_widget
            
            if self.manager.check_stage_status(stage_def.name):
                stage_widget.set_status(StageStatus.SUCCESS)
                artifacts = self.manager.get_stage_artifacts(stage_def.name)
                stage_widget.set_artifacts(artifacts)
            else:
                stage_widget.set_status(StageStatus.PENDING)

        self.run_all_btn.setEnabled(True)
        self.log(f"Pipeline initialized for {rcsb_id}. Loaded existing artifact status.")
        self.update_context_display()

    def run_stage(self, config: StageDefinition):
        if not self.manager:
            QMessageBox.warning(self, "Error", "Please initialize pipeline first")
            return
            
        if self.current_executor:
            QMessageBox.warning(self, "Busy", "Another stage is already running")
            return
            
        stage_widget = self.stage_widgets[config.name]
        params = stage_widget.get_params()
        
        self.log(f"Starting stage: {config.name} with params: {params}")
        stage_widget.set_status(StageStatus.RUNNING)
        self.update_context_display()
        
        self.current_executor = StageExecutor(self.manager, config.name, params)
        self.current_executor.finished.connect(self.on_stage_finished)
        self.current_executor.start()

    def on_stage_finished(self, stage_name: str, success: bool, result: Any, message: str):
        if stage_name not in self.stage_widgets:
            self.log(f"Received signal for unknown stage: {stage_name}", "ERROR")
            return

        stage_widget = self.stage_widgets[stage_name]

        if success:
            self.log(f"Stage completed: {stage_name}")
            stage_widget.set_status(StageStatus.SUCCESS)
            
            artifacts = self.manager.get_stage_artifacts(stage_name)
            stage_widget.set_artifacts(artifacts)
            
            self.view_stage_artifacts(stage_widget.config)
            
        else:
            self.log(f"Stage FAILED: {stage_name}\n{message}", "ERROR")
            stage_widget.set_status(StageStatus.FAILED)
            QMessageBox.critical(self, "Stage Failed", f"Stage {stage_name} failed:\n\n{message}")
            
        self.current_executor = None
        self.update_context_display()
        
        if hasattr(self, 'sequential_queue') and self.sequential_queue:
            if success:
                QTimer.singleShot(100, self.run_next_in_queue)
            else:
                self.log("Sequential run stopped due to failure.", "ERROR")
                self.sequential_queue.clear()

    # --- THIS IS THE ONLY METHOD MODIFIED in NPETPipelineOrchestrator ---
    def view_stage_artifacts(self, config: StageDefinition):
        """
        Called when 'View Artifacts' is clicked. 
        Gathers all necessary data and passes it to the viewer.
        """
        if config.name not in self.stage_widgets:
            self.log(f"Cannot view artifacts for unknown stage: {config.name}", "ERROR")
            return
        
        if not self.manager:
            self.log("Manager not initialized.", "ERROR")
            return

        stage_widget = self.stage_widgets[config.name]
        
        # Get the *current* parameters from the GUI
        current_params = stage_widget.get_params()
        
        # Pass the config, manager, and current params to the viewer
        self.artifact_viewer.set_stage_context(
            stage_config=config,
            manager=self.manager,
            current_params=current_params
        )
        # Switch to the artifact tab
        self.tabs.setCurrentWidget(self.artifact_viewer)
        
    def run_all_stages(self):
        if not self.manager:
            QMessageBox.warning(self, "Error", "Please initialize pipeline first")
            return
            
        self.sequential_queue = [stage.name for stage in self.manager.stages]
        self.run_next_in_queue()

    def run_next_in_queue(self):
        if not self.sequential_queue:
            self.log("All stages in queue completed.")
            return
            
        stage_name = self.sequential_queue.pop(0)
        stage_widget = self.stage_widgets[stage_name]
        
        if stage_widget.status == StageStatus.SUCCESS:
            self.log(f"Skipping already completed stage: {stage_name}")
            self.run_next_in_queue()
            return
            
        self.log(f"Running next in queue: {stage_name}")
        self.run_stage(stage_widget.config)
        
    def reset_pipeline(self):
        if self.current_executor:
            self.log("Cannot reset, a stage is running.", "WARN")
            return
            
        self.rcsb_id = None
        self.manager = None
        self.id_input.clear()
        self.run_all_btn.setEnabled(False)
        
        while self.stages_layout.count() > 1:
            child = self.stages_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.stage_widgets.clear()
        
        self.log_text.clear()
        self.context_table.setRowCount(0)
        self.artifact_viewer.plotter.clear()
        self.log("Pipeline reset. Enter an RCSB ID.")

    def update_context_display(self):
        self.context_table.setRowCount(0)
        if not self.manager:
            return
            
        self.context_table.setRowCount(len(self.manager.context))
        
        for i, (key, value) in enumerate(self.manager.context.items()):
            self.context_table.setItem(i, 0, QTableWidgetItem(str(key)))
            
            value_str = str(value)
            if isinstance(value, np.ndarray):
                value_str = f"np.ndarray(shape={value.shape}, dtype={value.dtype})"
            elif len(value_str) > 100:
                value_str = value_str[:100] + "..."
                
            self.context_table.setItem(i, 1, QTableWidgetItem(value_str))

def main():
    os.environ['QT_API'] = 'pyqt5'
    
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setStyle('Fusion') 
    
    orchestrator = NPETPipelineOrchestrator()
    orchestrator.show()
    
    sys.exit(app.exec_())

if __name__ == '__main__':
    if "RIBETL_DATA" not in os.environ:
        print("ERROR: 'RIBETL_DATA' environment variable not set.")
        # You should probably exit here, but we'll let it try
    
    pv.set_plot_theme("paraview") 
    pv.global_theme.multi_rendering_splitting_position = 0.5

    main()
```

pipeline_manager.py:
```
# Create a new file: ribctl/lib/npet/pipeline_manager.py

import os
import sys
from typing import TYPE_CHECKING

sys.path.append("/home/rtviii/dev/riboxyz")
import json
import numpy as np
import pyvista as pv
from pathlib import Path
from dataclasses import dataclass, field
from typing import Any, Dict, List, Callable, Optional

import open3d as o3d
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
import visualization_library as viz_lib

# Import your core computation functions
from ribctl.lib.npet.kdtree_approach import (
    landmark_constriction_site,
    ribosome_entities,
    filter_residues_parallel,
    transform_points_to_C0,
    create_point_cloud_mask,
    transform_points_from_C0,
    clip_pcd_via_ashape,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    apply_poisson_reconstruction,
)
from ribctl.lib.npet.alphalib import (
    cif_to_point_cloud,
    quick_surface_points,
    fast_normal_estimation,
    validate_mesh_pyvista,
)
from ribctl import RIBETL_DATA

if not RIBETL_DATA:
    raise EnvironmentError(
        "RIBETL_DATA environment variable is not set or not imported correctly."
    )

RIBETL_DATA_PATH = Path(RIBETL_DATA)


@dataclass
class VisualizationSpec:

    name: str

    function: Callable

    artifact_map: Dict[str, str] = field(default_factory=dict)
    context_map: Dict[str, str] = field(default_factory=dict)
    param_map: Dict[str, str] = field(default_factory=dict)


@dataclass
class StageDefinition:

    name: str
    description: str
    run_method_name: str  
    default_params: Dict[str, Any] = field(default_factory=dict)
    get_artifact_paths: Callable[[Path, Path], Dict[str, Path]] = field(
        default_factory=lambda: lambda a, b: {}
    )

    visualizations: List[VisualizationSpec] = field(default_factory=list)


class PipelineManager:
    """
    Orchestrates the NPET pipeline logic, decoupled from any UI.
    Manages state via the file system ("Artifact-First").
    """

    def __init__(self, rcsb_id: str):
        self.rcsb_id = rcsb_id.upper()
        self.assets_dir = RIBETL_DATA_PATH / self.rcsb_id
        self.artifacts_dir = self.assets_dir / "artifacts"
        self.artifacts_dir.mkdir(exist_ok=True, parents=True)

        self.stages = self._define_stages()
        self.stage_map = {s.name: s for s in self.stages}

        # This context holds *in-memory* data passed between stages *during a run*
        self.context: Dict[str, Any] = {"rcsb_id": self.rcsb_id}

    def get_stage_artifacts(self, stage_name: str) -> Dict[str, Path]:
        """Get the defined artifact paths for a given stage."""
        stage = self.stage_map[stage_name]
        return stage.get_artifact_paths(self.assets_dir, self.artifacts_dir)

    def check_stage_status(self, stage_name: str) -> bool:
        """Check if a stage is 'complete' by verifying its artifacts exist."""
        paths = self.get_stage_artifacts(stage_name).values()
        if not paths:
            return True  # Stages without artifacts (like Setup) are always "complete"
        return all(p.exists() for p in paths)

    def run_stage(self, stage_name: str, params: Dict[str, Any]) -> Dict[str, Any]:
        """
        Run a single stage by name.
        This is the function our QThread will call.
        """
        stage_def = self.stage_map[stage_name]

        # 1. Ensure dependencies (context) are met
        self._load_context_for_stage(stage_def)

        # 2. Run the actual compute method
        method = getattr(self, stage_def.run_method_name)
        stage_params = {**stage_def.default_params, **params}

        print(f"Running stage '{stage_name}' with params: {stage_params}")
        result = method(stage_params)

        # 4. Update the in-memory context
        if result:
            self.context.update(result)

        print(f"Stage '{stage_name}' complete.")
        return result

    def _load_context_for_stage(self, stage: StageDefinition):
        """
        Super simple dependency loading.
        If a context key is missing, try to load it from a predecessor's artifact.
        """
        # --- Dependencies for "PTC Identification" ---
        if "cifpath" not in self.context:
            setup_artifacts = self.get_stage_artifacts("Setup")
            cif_path = setup_artifacts.get("cif")
            if cif_path and cif_path.exists():
                self.context["cifpath"] = cif_path
                self.context["ro"] = RibosomeOps(self.rcsb_id)  # Re-create RibosomeOps
            else:
                if stage.name != "Setup":
                    raise FileNotFoundError(f"CIF file not found. Run 'Setup' first.")

        # --- Dependencies for "Entity Filtering" (and many others) ---
        if "ptc_pt" not in self.context and stage.name not in [
            "Setup",
            "PTC Identification",
        ]:
            ptc_artifact = self.get_stage_artifacts("PTC Identification").get(
                "ptc_json"
            )
            if ptc_artifact and ptc_artifact.exists():
                with open(ptc_artifact, "r") as f:
                    self.context["ptc_pt"] = np.array(json.load(f)["location"])
            else:
                raise FileNotFoundError(
                    f"PTC JSON not found. Run 'PTC Identification' first."
                )

        if "constriction_pt" not in self.context and stage.name not in [
            "Setup",
            "PTC Identification",
            "Constriction Site",
        ]:
            # This one is tricky, it's not a file. It's fetched.
            # Let's just fetch it if missing.
            try:
                self.context["constriction_pt"] = landmark_constriction_site(
                    self.rcsb_id
                )
            except Exception as e:
                raise RuntimeError(f"Could not load constriction point: {e}")

        # --- Dependencies for "Clustering" ---
        if "interior_points" not in self.context and stage.name in [
            "Clustering",
            "Refinement",
            "Surface Extraction",
        ]:
            pcd_artifact = self.get_stage_artifacts("Point Cloud Processing").get(
                "interior_points"
            )
            if pcd_artifact and pcd_artifact.exists():
                self.context["interior_points"] = np.load(pcd_artifact)
            else:
                raise FileNotFoundError(
                    f"Interior points not found. Run 'Point Cloud Processing' first."
                )

        # ... Add more context loaders as needed ...

    # ---------------------------------------------------------------------
    # STAGE DEFINITIONS
    # ---------------------------------------------------------------------

    def _define_stages(self) -> List[StageDefinition]:
        """
        This is the new "single source of truth" for the pipeline.
        We define all stages, their params, their artifacts, and their logic.
        """

        # This is a bit of a hack to avoid circular dependencies
        # A better design would be to put viz funcs in a new `visualization_library.py`
        # and have both the manager and orchestrator import from it.
        # But for now, this fulfills the "all in one file" request.
        try:
            from npet_orchestrator import (
                viz_simple_mesh,
                viz_landmarks_and_pointcloud,
                viz_single_landmark,
            )
        except ImportError:
            print(
                "WARNING [PipelineManager]: Could not import visualization functions from npet_orchestrator."
            )
            print("This is OK during initialization, but will fail if GUI is run.")

            # Define dummy functions to allow the file to load
            def viz_simple_mesh(*args, **kwargs):
                pass

            def viz_landmarks_and_pointcloud(*args, **kwargs):
                pass

            def viz_single_landmark(*args, **kwargs):
                pass

        return [
            StageDefinition(
                name="Setup",
                description="Initialize pipeline and load structure data",
                run_method_name="_run_setup",
                get_artifact_paths=lambda assets, artifacts: {
                    "cif": assets / f"{self.rcsb_id}.cif"
                },
                visualizations=[],  # No 3D viz for this stage
            ),
            StageDefinition(
                name="PTC Identification",
                description="Identify Peptidyl Transferase Center location",
                run_method_name="_run_ptc_identification",
                get_artifact_paths=lambda assets, artifacts: {
                    "ptc_json": assets / f"{self.rcsb_id}_PTC.json"
                },
                visualizations=[
                    VisualizationSpec(
                        name="View PTC Landmark",
                        function=viz_lib.viz_single_landmark,  # <-- Use viz_lib
                        artifact_map={"landmark_json": "ptc_json"},
                        context_map={"rcsb_id": "rcsb_id"},
                    )
                ],
            ),
            StageDefinition(
                name="Constriction Site",
                description="Identify ribosome exit tunnel constriction site (fetches from API)",
                run_method_name="_run_constriction_identification",
                visualizations=[
                    VisualizationSpec(
                        name="View PTC + Constriction",
                        function=viz_lib.viz_landmarks_and_pointcloud,  # <-- Use viz_lib
                        artifact_map={"background_mesh": "ashape_mesh"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Alpha Shape (Exterior)",
                description="Generate exterior ribosome surface mesh",
                run_method_name="_run_alpha_shape",
                default_params={
                    "d3d_alpha": 200,
                    "d3d_tol": 10,
                    "d3d_offset": 3,
                    "kdtree_radius": 40,
                    "max_nn": 60,
                    "tangent_planes_k": 20,
                    "PR_depth": 6,
                    "PR_ptweight": 4,
                },
                get_artifact_paths=lambda assets, artifacts: {
                    "ashape_mesh": assets / f"{self.rcsb_id}_ALPHA_SHAPE.ply",
                    "ashape_mesh_ascii": assets
                    / f"{self.rcsb_id}_ALPHA_SHAPE_ascii.ply",
                    "structure_ptcloud": artifacts
                    / f"{self.rcsb_id}_structure_ptcloud.npy",
                    "surface_points": artifacts
                    / f"{self.rcsb_id}_alpha_surface_points.npy",
                    "normals": artifacts / f"{self.rcsb_id}_alpha_normals.ply",
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Exterior Mesh",
                        function=viz_simple_mesh,
                        artifact_map={"mesh_to_show": "ashape_mesh"},
                        context_map={"rcsb_id": "rcsb_id"},
                    )
                ],
            ),
            StageDefinition(
                name="Entity Filtering",
                description="Filter atoms within tunnel cylinder region",
                run_method_name="_run_entity_filtering",
                default_params={"radius": 35, "height": 120},
                get_artifact_paths=lambda assets, artifacts: {
                    "filtered_points": artifacts / f"{self.rcsb_id}_filtered_points.npy"
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Filtered Points + Landmarks",
                        function=viz_landmarks_and_pointcloud,
                        artifact_map={"point_cloud_data": "filtered_points"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Point Cloud Processing",
                description="Transform points, create mask, get interior",
                run_method_name="_run_point_cloud_processing",
                default_params={"voxel_size": 1, "atom_size": 2},
                get_artifact_paths=lambda assets, artifacts: {
                    "transformed_points": artifacts
                    / f"{self.rcsb_id}_transformed_points.npy",
                    "interior_points": artifacts
                    / f"{self.rcsb_id}_interior_points.npy",
                    "back_projected_points": artifacts
                    / f"{self.rcsb_id}_back_projected.npy",
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Interior Points + Landmarks",
                        function=viz_landmarks_and_pointcloud,
                        artifact_map={"point_cloud_data": "interior_points"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Clustering",
                description="Cluster tunnel points using DBSCAN",
                run_method_name="_run_clustering",
                default_params={"epsilon": 5.5, "min_samples": 600},
                get_artifact_paths=lambda assets, artifacts: {
                    "largest_cluster": artifacts
                    / f"{self.rcsb_id}_largest_cluster.npy",
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Largest Cluster",
                        function=viz_landmarks_and_pointcloud,
                        artifact_map={"point_cloud_data": "largest_cluster"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Refinement",
                description="Refine clusters with secondary DBSCAN",
                run_method_name="_run_refinement",
                default_params={"epsilon": 3.5, "min_samples": 175},
                get_artifact_paths=lambda assets, artifacts: {
                    "refined_cluster": artifacts / f"{self.rcsb_id}_refined_cluster.npy"
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Refined Cluster",
                        function=viz_landmarks_and_pointcloud,
                        artifact_map={"point_cloud_data": "refined_cluster"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Surface Extraction",
                description="Extract surface points from refined cluster",
                run_method_name="_run_surface_extraction",
                default_params={"alpha": 2, "tolerance": 1, "offset": 2},
                get_artifact_paths=lambda assets, artifacts: {
                    "surface_points": artifacts / f"{self.rcsb_id}_surface_points.npy"
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Surface Points",
                        function=viz_landmarks_and_pointcloud,
                        artifact_map={"point_cloud_data": "surface_points"},
                        context_map={
                            "ptc_coord": "ptc_pt",
                            "constriction_coord": "constriction_pt",
                            "rcsb_id": "rcsb_id",
                        },
                    )
                ],
            ),
            StageDefinition(
                name="Normal Estimation",
                description="Estimate normals for surface reconstruction",
                run_method_name="_run_normal_estimation",
                default_params={
                    "kdtree_radius": 10,
                    "kdtree_max_nn": 15,
                    "correction_tangent_planes_n": 10,
                },
                get_artifact_paths=lambda assets, artifacts: {
                    "normal_estimated_pcd": artifacts
                    / f"{self.rcsb_id}_normal_estimated_pcd.ply"
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Normals PCD",
                        function=viz_simple_mesh,
                        artifact_map={"mesh_to_show": "normal_estimated_pcd"},
                        context_map={"rcsb_id": "rcsb_id"},
                    )
                ],
            ),
            StageDefinition(
                name="Mesh Reconstruction",
                description="Reconstruct tunnel mesh using Poisson",
                run_method_name="_run_mesh_reconstruction",
                default_params={"depth": 6, "ptweight": 3},
                get_artifact_paths=lambda assets, artifacts: {
                    "npet_mesh": assets / f"{self.rcsb_id}_NPET_MESH.ply",
                    "npet_mesh_ascii": assets / f"{self.rcsb_id}_NPET_MESH_ascii.ply",
                },
                visualizations=[
                    VisualizationSpec(
                        name="View Final NPET Mesh",
                        function=viz_simple_mesh,
                        artifact_map={"mesh_to_show": "npet_mesh"},
                        context_map={"rcsb_id": "rcsb_id"},
                    )
                ],
            ),
        ]

    # ---------------------------------------------------------------------
    # --- STAGE COMPUTE METHODS (Your existing logic) ---
    # ---------------------------------------------------------------------

    def _run_setup(self, params: Dict) -> Dict[str, Any]:
        """Ported from SetupStage"""

        cif_path = self.get_stage_artifacts("Setup").get("cif")
        if not cif_path:
            raise FileNotFoundError(f"CIF path object is None! Check StageDefinition.")

        if not cif_path.exists():
            parent_dir = cif_path.parent
            print(
                f"DEBUG [_run_setup]: Parent dir '{parent_dir}' exists: {parent_dir.exists()}"
            )
            if parent_dir.exists():
                print(
                    f"DEBUG [_run_setup]: Contents of {parent_dir}: {os.listdir(parent_dir)}"
                )
            raise FileNotFoundError(f"CIF file not found at {cif_path}")

        ro = RibosomeOps(self.rcsb_id)

        return {
            "cifpath": cif_path,
            "ro": ro,
        }

    def _run_ptc_identification(self, params: Dict) -> Dict[str, Any]:
        """Ported from PTCIdentificationStage"""
        ptc_info = PTC_location(self.rcsb_id)
        ptc_pt = np.array(ptc_info.location)

        ptc_json_path = self.get_stage_artifacts("PTC Identification")["ptc_json"]
        with open(ptc_json_path, "w") as f:
            f.write(ptc_info.model_dump_json(indent=2))

        return {"ptc_info": ptc_info, "ptc_pt": ptc_pt}

    def _run_constriction_identification(self, params: Dict) -> Dict[str, Any]:
        """Ported from ConstrictionIdentificationStage"""
        try:
            constriction_pt = landmark_constriction_site(self.rcsb_id)
            return {"constriction_pt": constriction_pt}
        except Exception as e:
            raise RuntimeError(f"Failed to fetch constriction site: {e}")

    def _run_alpha_shape(self, params: Dict) -> Dict[str, Any]:
        """Ported from alphalib.py"""
        cif_path = self.context["cifpath"]
        ro = self.context["ro"]

        ptcloud_path = self.get_stage_artifacts("Alpha Shape (Exterior)")[
            "structure_ptcloud"
        ]
        surface_pts_path = self.get_stage_artifacts("Alpha Shape (Exterior)")[
            "surface_points"
        ]
        normals_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["normals"]
        mesh_path = self.get_stage_artifacts("Alpha Shape (Exterior)")["ashape_mesh"]

        print("Extracting point cloud from CIF...")
        chains = ro.first_assembly_auth_asym_ids()
        ptcloud = cif_to_point_cloud(cif_path, chains, do_atoms=True)
        np.save(ptcloud_path, ptcloud)

        print("Calculating surface points (Delaunay 3D)...")
        surface_pts = quick_surface_points(
            ptcloud, params["d3d_alpha"], params["d3d_tol"], params["d3d_offset"]
        )
        np.save(surface_pts_path, surface_pts)

        print("Estimating normals...")
        normal_estimated_pcd = fast_normal_estimation(
            surface_pts,
            params["kdtree_radius"],
            params["max_nn"],
            params["tangent_planes_k"],
        )
        o3d.io.write_point_cloud(str(normals_path), normal_estimated_pcd)

        print("Applying Poisson Reconstruction...")
        apply_poisson_reconstruction(
            str(normals_path),
            mesh_path,
            recon_depth=params["PR_depth"],
            recon_pt_weight=params["PR_ptweight"],
        )

        # Post-processing: keep largest component
        mesh = pv.read(mesh_path)
        labeled = mesh.connectivity(largest=True)
        labeled.save(mesh_path)

        watertight = validate_mesh_pyvista(labeled)
        print(f"Alpha shape mesh generated. Watertight: {watertight}")

        return {"ashape_mesh": mesh_path, "watertight": watertight}

    def _run_entity_filtering(self, params: Dict) -> Dict[str, Any]:
        """Ported from EntityFilteringStage"""
        cifpath = self.context["cifpath"]
        ptc_pt = self.context["ptc_pt"]
        constriction_pt = self.context["constriction_pt"]

        # Logic for tunnel_debris (hardcoded)
        tunnel_debris = {
            "3J7Z": ["a", "7"],
            "5GAK": ["z"],
            "5NWY": ["s"],
            "7A5G": ["Y2"],
            "9F1D": ["BK"],
        }
        skip_chains = tunnel_debris.get(self.rcsb_id, [])

        residues = ribosome_entities(self.rcsb_id, cifpath, "R", skip_chains)

        filtered_residues = filter_residues_parallel(
            residues, ptc_pt, constriction_pt, params["radius"], params["height"]
        )

        filtered_points = np.array(
            [
                atom.get_coord()
                for residue in filtered_residues
                for atom in residue.child_list
            ]
        )

        filtered_points_path = self.get_stage_artifacts("Entity Filtering")[
            "filtered_points"
        ]
        np.save(filtered_points_path, filtered_points)

        return {
            "filtered_residues": filtered_residues,
            "filtered_points": filtered_points,
        }

    def _run_point_cloud_processing(self, params: Dict) -> Dict[str, Any]:
        """Ported from PointCloudProcessingStage"""
        filtered_points = self.context["filtered_points"]
        ptc_pt = self.context["ptc_pt"]
        constriction_pt = self.context["constriction_pt"]
        ashape_mesh_path = self.get_stage_artifacts("Alpha Shape (Exterior)")[
            "ashape_mesh"
        ]

        # 1. Transform points to cylinder coordinate system
        transformed_points = transform_points_to_C0(
            filtered_points, ptc_pt, constriction_pt
        )

        # 2. Create voxel mask
        mask, (x, y, z) = create_point_cloud_mask(
            transformed_points,
            radius=self.stage_map["Entity Filtering"].default_params["radius"],
            height=self.stage_map["Entity Filtering"].default_params["height"],
            voxel_size=params["voxel_size"],
            radius_around_point=params["atom_size"],
        )

        # 3. Get "empty" voxel centers
        X, Y, Z = np.meshgrid(x, y, z, indexing="ij")
        empty_voxels_C0 = np.column_stack((X[~mask], Y[~mask], Z[~mask]))

        # 4. Transform "empty" points back to world coordinates
        empty_in_world_coords = transform_points_from_C0(
            empty_voxels_C0, ptc_pt, constriction_pt
        )

        # 5. Clip empty points by the exterior alpha shape
        ashape_mesh = pv.read(ashape_mesh_path)
        interior_points, _ = clip_pcd_via_ashape(empty_in_world_coords, ashape_mesh)

        # Save artifacts
        np.save(
            self.get_stage_artifacts("Point Cloud Processing")["transformed_points"],
            transformed_points,
        )
        np.save(
            self.get_stage_artifacts("Point Cloud Processing")["interior_points"],
            interior_points,
        )
        np.save(
            self.get_stage_artifacts("Point Cloud Processing")["back_projected_points"],
            empty_in_world_coords,
        )

        return {"interior_points": interior_points}

    def _run_clustering(self, params: Dict) -> Dict[str, Any]:
        """Ported from ClusteringStage"""
        interior_points = self.context["interior_points"]

        db, clusters_container = DBSCAN_capture(
            interior_points, params["epsilon"], params["min_samples"]
        )

        largest_cluster_pts, cluster_id = DBSCAN_pick_largest_cluster(
            clusters_container
        )

        print(
            f"Picked largest cluster #{cluster_id} with {len(largest_cluster_pts)} points."
        )

        # Save artifacts
        np.save(
            self.get_stage_artifacts("Clustering")["largest_cluster"],
            largest_cluster_pts,
        )

        return {"largest_cluster": largest_cluster_pts}

    def _run_refinement(self, params: Dict) -> Dict[str, Any]:
        """Ported from RefinementStage"""
        largest_cluster = self.context["largest_cluster"]

        db, clusters_container = DBSCAN_capture(
            largest_cluster, params["epsilon"], params["min_samples"]
        )

        refined_cluster_pts, cluster_id = DBSCAN_pick_largest_cluster(
            clusters_container
        )

        print(
            f"Picked largest refined cluster #{cluster_id} with {len(refined_cluster_pts)} points."
        )

        np.save(
            self.get_stage_artifacts("Refinement")["refined_cluster"],
            refined_cluster_pts,
        )

        return {"refined_cluster": refined_cluster_pts}

    def _run_surface_extraction(self, params: Dict) -> Dict[str, Any]:
        """Ported from SurfaceExtractionStage"""
        refined_cluster = self.context["refined_cluster"]

        surface_points = quick_surface_points(
            refined_cluster, params["alpha"], params["tolerance"], params["offset"]
        )

        np.save(
            self.get_stage_artifacts("Surface Extraction")["surface_points"],
            surface_points,
        )

        return {"surface_points": surface_points}

    def _run_normal_estimation(self, params: Dict) -> Dict[str, Any]:
        """Ported from NormalEstimationStage"""
        surface_points = self.context["surface_points"]

        normal_estimated_pcd = fast_normal_estimation(
            surface_points,
            params["kdtree_radius"],
            params["kdtree_max_nn"],
            params["correction_tangent_planes_n"],
        )

        pcd_path = self.get_stage_artifacts("Normal Estimation")["normal_estimated_pcd"]
        o3d.io.write_point_cloud(str(pcd_path), normal_estimated_pcd)

        return {"normal_estimated_pcd_path": pcd_path}

    def _run_mesh_reconstruction(self, params: Dict) -> Dict[str, Any]:
        """Ported from MeshReconstructionStage"""
        pcd_path = self.context["normal_estimated_pcd_path"]
        mesh_path = self.get_stage_artifacts("Mesh Reconstruction")["npet_mesh"]

        apply_poisson_reconstruction(
            str(pcd_path),
            mesh_path,
            recon_depth=params["depth"],
            recon_pt_weight=params["ptweight"],
        )

        mesh = pv.read(mesh_path)
        labeled = mesh.connectivity(largest=True)
        labeled.save(mesh_path)

        watertight = validate_mesh_pyvista(labeled)
        print(f"NPET mesh generated. Watertight: {watertight}")

        return {"npet_mesh": mesh_path, "npet_watertight": watertight}
```


visualization_library.py:
```
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
```

Tell me if you need to see any other files and also if you want me to disambiguate any specification from what i just described.