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
from typing import Dict, Any, Optional, List
from enum import Enum
import pickle
from datetime import datetime

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
from PyQt5.QtGui import QFont, QColor, QPalette

from pipeline_manager import PipelineManager, StageDefinition

class StageStatus(Enum):
    PENDING = "pending"
    RUNNING = "running"
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped" # We can keep this for manual override

class StageExecutor(QThread):
    """
    Worker thread for running a pipeline stage.
    This fixes the "jank" by running compute off the main thread.
    """
    # Signal: stage_name, success, result_dict, error_message
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

# --- StageWidget (Modified) ---
class StageWidget(QFrame):
    """Widget representing a single pipeline stage"""
    
    run_requested = pyqtSignal(StageDefinition)
    view_artifacts = pyqtSignal(StageDefinition)
    
    def __init__(self, config: StageDefinition, parent=None): # <-- Takes new StageDefinition
        super().__init__(parent)
        self.config = config
        self.status = StageStatus.PENDING
        self.params = config.default_params.copy()
        
        # This now holds the *paths* to the artifacts
        self.artifacts: Dict[str, Path] = {} 
        
        self.execution_time = None
        self.setup_ui()
        self.update_status_display()
        
    def setup_ui(self):
        # ... (Your existing setup_ui is fine) ...
        # (Make sure to replace config.stage_class with config.name where appropriate)
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
        # ... (This is mostly the same) ...
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
        # ... (This is the same) ...
        params = {}
        for name, widget in self.param_widgets.items():
            if isinstance(widget, QComboBox):
                params[name] = widget.currentText() == "True"
            elif isinstance(widget, QSpinBox):
                params[name] = widget.value()
            elif isinstance(widget, QDoubleSpinBox):
                params[name] = widget.value()
            else: # QLineEdit
                # Attempt to auto-cast str to float/int
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
        """Store the artifact paths and update the UI."""
        self.artifacts = artifacts
        self.update_status_display()

# --- ArtifactViewer (Modified) ---
class ArtifactViewer(QWidget):
    
    def __init__(self, parent=None):
        super().__init__(parent)
        self.current_artifacts: Dict[str, Path] = {}
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout(self)
        
        info_group = QGroupBox("Artifact Info")
        info_layout = QVBoxLayout()
        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setMaximumHeight(150)
        self.info_text.setFont(QFont("Courier New", 9))
        info_layout.addWidget(self.info_text)
        info_group.setLayout(info_layout)
        layout.addWidget(info_group)
        
        selector_layout = QHBoxLayout()
        selector_layout.addWidget(QLabel("Artifact:"))
        self.artifact_combo = QComboBox()
        self.artifact_combo.currentIndexChanged.connect(self.load_selected_artifact)
        selector_layout.addWidget(self.artifact_combo, 1)
        self.export_btn = QPushButton("Export...")
        self.export_btn.clicked.connect(self.export_artifact)
        selector_layout.addWidget(self.export_btn)
        layout.addLayout(selector_layout)
        
        # --- REVERT TO THIS ---
        # 3D viewer
        self.plotter_frame = QFrame()
        self.plotter_frame.setFrameStyle(QFrame.StyledPanel | QFrame.Sunken)
        plotter_layout = QVBoxLayout(self.plotter_frame)
        
        # Create the plotter with the frame as its parent
        self.plotter = QtInteractor(self.plotter_frame)
        
        # Tell the plotter to use 3D controls
        self.plotter.enable_trackball_style()
        
        # Add the plotter to the frame's layout
        plotter_layout.addWidget(self.plotter)  
        
        # Add the frame (which contains the plotter) to the main layout
        layout.addWidget(self.plotter_frame, 1)
        # --- END REVERT ---
        
    def set_artifacts(self, artifacts: Dict[str, Path], stage_name: str):
        self.current_artifacts = artifacts
        self.artifact_combo.clear()
        
        for name, path in artifacts.items():
            if path.exists():
                # Add the logical name and store the path object as data
                self.artifact_combo.addItem(f"{name} ({path.name})", userData=path)
            
        self.info_text.setPlainText(f"Stage: {stage_name}\n{len(artifacts)} artifact(s) defined.")
        
        if self.artifact_combo.count() > 0:
            self.load_selected_artifact(0)
        else:
            self.plotter.clear()
            self.info_text.append("\nNo viewable artifacts found on disk.")

    def load_selected_artifact(self, index: int):
        if index < 0:
            self.plotter.clear()
            return
            
        artifact_path: Path = self.artifact_combo.itemData(index)
        if not artifact_path or not artifact_path.exists():
            self.info_text.append(f"\nERROR: File not found: {artifact_path}")
            return
            
        self.plotter.clear()
        
        # *** FIX 1: Typo fix ***
        # 'enable_zoom_scaling' was wrong, it's 'enable_zoom_style'
        self.plotter.enable_trackball_style()
        # *** END FIX 1 ***
        
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
            
            # *** FIX 2: Add handler for .cif files ***
            elif ext == '.cif':
                info = f"Structure File: {artifact_path.name}\n"
                info += f"Size: {artifact_path.stat().st_size / (1024*1024):.2f} MB\n\n"
                info += "Preview not available in this viewer."
                self.info_text.setPlainText(info)
                self.plotter.add_text("CIF structure file (see info panel)", position='upper_left') # <-- FIX
            # *** END FIX 2 ***

            elif ext == '.json':
                with open(artifact_path) as f:
                    data = json.load(f)
                self.info_text.setPlainText(json.dumps(data, indent=2))
                self.plotter.add_text("JSON data (see info panel)", position='upper_left') # <-- FIX
                
            else:
                self.info_text.setPlainText(f"Cannot preview file type: {ext}")
                self.plotter.add_text(f"Cannot preview {ext}", position='upper_left') # <-- FIX

            self.plotter.reset_camera()
            
        except Exception as e:
            self.info_text.setPlainText(f"Error loading artifact:\n{str(e)}")
            print(traceback.format_exc())

    def export_artifact(self):
        # ... (This is the same) ...
        pass

class NPETPipelineOrchestrator(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.rcsb_id: Optional[str] = None
        self.manager: Optional[PipelineManager] = None
        self.stage_widgets: Dict[str, StageWidget] = {}
        self.current_executor: Optional[StageExecutor] = None
        
        # This is no longer needed, it's in the manager
        # self.stages = self.define_stages() 
        
        self.setup_ui()
        self.create_menu_bar()
        
        self.setWindowTitle("NPET Pipeline Orchestrator")
        self.setGeometry(100, 100, 1400, 900)
        self.log("Orchestrator started. Enter an RCSB ID to begin.")

    def create_menu_bar(self):
        # ... (Your menubar code is fine, omitted for brevity) ...
        pass

    # ... (Your dialog functions: run_from_stage_dialog, etc. are fine) ...
    # ... (They will need to be adapted to use stage *names* or *indices* from self.manager.stages) ...

    def define_stages(self) -> List[StageDefinition]:
        """DEPRECATED: This logic is now in PipelineManager."""
        # This function is no longer called.
        # We get stage definitions from self.manager.stages
        return []
        
    def setup_ui(self):
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        
        # Top control bar
        control_layout = QHBoxLayout()
        control_layout.addWidget(QLabel("RCSB ID:"))
        self.id_input = QLineEdit()
        self.id_input.setPlaceholderText("e.g., 7K00")
        self.id_input.setMaximumWidth(100)
        self.id_input.returnPressed.connect(self.initialize_pipeline) # Hook up Enter key
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
        
        # Main splitter
        splitter = QSplitter(Qt.Horizontal)
        
        # Left: Stage list
        stages_widget = QWidget()
        stages_layout_v = QVBoxLayout(stages_widget)
        stages_label = QLabel("Pipeline Stages")
        stages_label.setFont(QFont("Arial", 14, QFont.Bold))
        stages_layout_v.addWidget(stages_label)
        
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setMinimumWidth(400)
        
        stages_container = QWidget()
        # This layout will be populated by initialize_pipeline
        self.stages_layout = QVBoxLayout(stages_container) 
        self.stages_layout.addStretch() # Add stretch at the bottom
        scroll.setWidget(stages_container)
        stages_layout_v.addWidget(scroll)
        splitter.addWidget(stages_widget)
        
        # Right: Tabs
        tabs = QTabWidget()
        self.artifact_viewer = ArtifactViewer()
        tabs.addTab(self.artifact_viewer, "Artifact Viewer")
        
        # ... (Your Log tab and Context tab setup is fine) ...
        log_widget = QWidget()
        log_layout = QVBoxLayout(log_widget)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setFont(QFont("Courier New", 9))
        log_layout.addWidget(self.log_text)
        tabs.addTab(log_widget, "Execution Log")
        
        context_widget = QWidget()
        context_layout = QVBoxLayout(context_widget)
        self.context_table = QTableWidget()
        self.context_table.setColumnCount(2)
        self.context_table.setHorizontalHeaderLabels(["Key", "Value"])
        self.context_table.horizontalHeader().setStretchLastSection(True)
        context_layout.addWidget(self.context_table)
        tabs.addTab(context_widget, "Pipeline Context")
        
        splitter.addWidget(tabs)
        splitter.setSizes([450, 950]) # Adjust initial split
        main_layout.addWidget(splitter)

    def log(self, message: str, level: str = "INFO"):
        timestamp = datetime.now().strftime("%H:%M:%S")
        formatted = f"[{timestamp}] {level}: {message}"
        self.log_text.append(formatted)
        print(formatted) # Also print to console
        
    def initialize_pipeline(self):
        """
        *** THIS IS THE KEY CHANGE ***
        Initializes the manager and loads state from the filesystem.
        """
        rcsb_id = self.id_input.text().strip().upper()
        if not rcsb_id:
            QMessageBox.warning(self, "Error", "Please enter an RCSB ID")
            return
            
        self.log(f"Initializing pipeline for {rcsb_id}...")
        self.rcsb_id = rcsb_id
        
        try:
            self.manager = PipelineManager(self.rcsb_id)
        except Exception as e:
            self.log(f"Failed to initialize manager: {e}", "ERROR")
            QMessageBox.critical(self, "Error", f"Failed to initialize manager:\n{e}")
            return

        # Clear old stage widgets
        while self.stages_layout.count() > 1: # Keep the stretch
            child = self.stages_layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()
        self.stage_widgets.clear()

        # Create new stage widgets from manager's definitions
        for stage_def in self.manager.stages:
            stage_widget = StageWidget(stage_def)
            stage_widget.run_requested.connect(self.run_stage)
            stage_widget.view_artifacts.connect(self.view_stage_artifacts)
            
            # Insert *before* the stretch
            self.stages_layout.insertWidget(self.stages_layout.count() - 1, stage_widget)
            self.stage_widgets[stage_def.name] = stage_widget
            
            # *** CHECK FILE SYSTEM FOR ARTIFACTS ***
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
        """
        *** THIS IS THE KEY CHANGE ***
        Runs a single stage using the QThread executor.
        """
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
        self.update_context_display() # Clear old context
        
        # Run in thread
        self.current_executor = StageExecutor(self.manager, config.name, params)
        self.current_executor.finished.connect(self.on_stage_finished)
        self.current_executor.start()

    def on_stage_finished(self, stage_name: str, success: bool, result: Any, message: str):
        """
        *** THIS IS THE KEY CHANGE ***
        Handles the 'finished' signal from the StageExecutor thread.
        """
        if stage_name not in self.stage_widgets:
            self.log(f"Received signal for unknown stage: {stage_name}", "ERROR")
            return

        stage_widget = self.stage_widgets[stage_name]

        if success:
            self.log(f"Stage completed: {stage_name}")
            stage_widget.set_status(StageStatus.SUCCESS)
            
            # Get artifact paths from manager and update widget
            artifacts = self.manager.get_stage_artifacts(stage_name)
            stage_widget.set_artifacts(artifacts)
            
            # Automatically view the new artifact
            self.view_stage_artifacts(stage_widget.config)
            
        else:
            self.log(f"Stage FAILED: {stage_name}\n{message}", "ERROR")
            stage_widget.set_status(StageStatus.FAILED)
            QMessageBox.critical(self, "Stage Failed", f"Stage {stage_name} failed:\n\n{message}")
            
        self.current_executor = None
        self.update_context_display()
        
        # Handle sequential run
        if hasattr(self, 'sequential_queue') and self.sequential_queue:
            if success:
                QTimer.singleShot(100, self.run_next_in_queue) # Run next on success
            else:
                self.log("Sequential run stopped due to failure.", "ERROR")
                self.sequential_queue.clear()

    def view_stage_artifacts(self, config: StageDefinition):
        stage_widget = self.stage_widgets[config.name]
        self.artifact_viewer.set_artifacts(stage_widget.artifacts, config.name)
        
    def run_all_stages(self):
        if not self.manager:
            QMessageBox.warning(self, "Error", "Please initialize pipeline first")
            return
            
        # Create a queue of stage names to run
        self.sequential_queue = [stage.name for stage in self.manager.stages]
        self.run_next_in_queue()

    def run_next_in_queue(self):
        if not self.sequential_queue:
            self.log("All stages in queue completed.")
            return
            
        stage_name = self.sequential_queue.pop(0)
        stage_widget = self.stage_widgets[stage_name]
        
        # Skip stages that are already successful
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
        
        # Clear stage widgets
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
    # CRITICAL: Set this BEFORE importing pyvista or creating QApplication
    # This is a good practice for PyVista/Qt compatibility
    os.environ['QT_API'] = 'pyqt5'
    
    # Enable high-DPI scaling
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setStyle('Fusion') # A clean, modern style
    
    orchestrator = NPETPipelineOrchestrator()
    orchestrator.show()
    
    sys.exit(app.exec_())

if __name__ == '__main__':
    # Make sure RIBETL_DATA is set
    if "RIBETL_DATA" not in os.environ:
        print("ERROR: 'RIBETL_DATA' environment variable not set.")
        print("Please set it to your main data directory.")
        # As a fallback for testing, we can set it here:
        # os.environ["RIBETL_DATA"] = "/path/to/your/RIBETL_DATA"
        # sys.exit(1)
        # For this example, I'll *assume* it's set, but this check is good practice.
    
    # This is also critical for PyVista/PyQt on some systems
    pv.set_plot_theme("paraview") # <-- Use this instead
    pv.global_theme.multi_rendering_splitting_position = 0.5

    main()