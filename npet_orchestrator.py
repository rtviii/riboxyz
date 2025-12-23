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
from dataclasses import dataclass, field
import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor
from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QTextEdit,
    QSplitter,
    QScrollArea,
    QFrame,
    QLineEdit,
    QGroupBox,
    QGridLayout,
    QSpinBox,
    QDoubleSpinBox,
    QComboBox,
    QFileDialog,
    QProgressBar,
    QMessageBox,
    QTabWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer
from PyQt5.QtGui import QFont

from pipeline_manager import PipelineManager, StageDefinition, VisualizationSpec


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

    def __init__(
        self, manager: PipelineManager, stage_name: str, params: Dict[str, Any]
    ):
        super().__init__()
        self.manager = manager
        self.stage_name = stage_name
        self.params = params

    def run(self):
        try:
            result = self.manager.run_stage(self.stage_name, self.params)
            self.finished.emit(
                self.stage_name, True, result, "Stage completed successfully"
            )
        except Exception as e:
            error_msg = f"Stage failed: {str(e)}\n{traceback.format_exc()}"
            self.finished.emit(self.stage_name, False, None, error_msg)


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
        """Simple method to set artifacts and update display."""
        self.artifacts = artifacts
        self.update_status_display()

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

        # NEW: Enable the ribosome toggle immediately
        self.artifact_viewer.set_manager(self.manager)

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

    def on_stage_finished(
        self, stage_name: str, success: bool, result: Any, message: str
    ):
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
            QMessageBox.critical(
                self, "Stage Failed", f"Stage {stage_name} failed:\n\n{message}"
            )

        self.current_executor = None
        self.update_context_display()

        if hasattr(self, "sequential_queue") and self.sequential_queue:
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
            stage_config=config, manager=self.manager, current_params=current_params
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


class ArtifactViewer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.manager: Optional[PipelineManager] = None
        self.stage_config: Optional[StageDefinition] = None
        self.current_params: Dict[str, Any] = {}
        self.all_artifacts: Dict[str, Path] = {}
        self.show_ribosome_structure = False  # NEW: Toggle state

        self.setup_ui()

    def set_manager(self, manager: PipelineManager):
        """Set the manager and enable global controls."""
        self.manager = manager
        self.ribosome_toggle.setEnabled(True)
        
        # Pre-load context for Setup stage to get CIF path
        try:
            setup_stage = next(s for s in manager.stages if s.name == "Setup")
            manager._load_context_for_stage(setup_stage)
        except Exception as e:
            print(f"Warning: Could not pre-load Setup context: {e}")

    def setup_ui(self):
        layout = QVBoxLayout(self)

        # NEW: Add global toggles section at the top
        toggles_group = QGroupBox("Global Overlays")
        toggles_layout = QHBoxLayout()

        self.ribosome_toggle = QPushButton("Show Ribosome Structure")
        self.ribosome_toggle.setCheckable(True)
        self.ribosome_toggle.setChecked(False)
        self.ribosome_toggle.clicked.connect(self.toggle_ribosome_structure)
        self.ribosome_toggle.setEnabled(False)  # Enable after stage is set
        toggles_layout.addWidget(self.ribosome_toggle)

        # Placeholder for future toggles
        toggles_layout.addStretch()
        toggles_group.setLayout(toggles_layout)
        layout.addWidget(toggles_group)

        # Info group
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

    def toggle_ribosome_structure(self):
        """NEW: Toggle method for ribosome structure"""
        self.show_ribosome_structure = self.ribosome_toggle.isChecked()

        if self.show_ribosome_structure:
            self.add_ribosome_overlay()
        else:
            self.remove_ribosome_overlay()

    def add_ribosome_overlay(self):
        """NEW: Add ribosome overlay"""
        if not self.manager:
            return

        try:
            # Import the viz function
            from visualization_library import viz_ribosome_structure

            # Get CIF path from context or artifacts
            cif_path = self.manager.context.get("cifpath")
            if not cif_path:
                setup_artifacts = self.manager.get_stage_artifacts("Setup")
                cif_path = setup_artifacts.get("cif")

            if not cif_path or not Path(cif_path).exists():
                self.info_text.append("\nError: CIF file not found")
                self.ribosome_toggle.setChecked(False)
                self.show_ribosome_structure = False
                return

            # Get chains from RibosomeOps
            ro = self.manager.context.get("ro")
            chains = ro.first_assembly_auth_asym_ids() if ro else None

            viz_ribosome_structure(
                self.plotter,
                str(cif_path),
                self.manager.rcsb_id,
                chains_to_show=chains,
                alpha=0.2,
                point_size=2,
            )

            self.info_text.append("\nRibosome structure overlay: ON")

        except Exception as e:
            self.info_text.append(f"\nFailed to add ribosome overlay: {e}")
            self.ribosome_toggle.setChecked(False)
            self.show_ribosome_structure = False

    def remove_ribosome_overlay(self):
        """NEW: Remove ribosome overlay"""
        try:
            # Remove by actor name
            self.plotter.remove_actor("ribosome_structure")
            self.plotter.remove_actor("ribosome_label")
            self.info_text.append("\nRibosome structure overlay: OFF")
        except:
            pass  # Actor might not exist

    def set_stage_context(
        self,
        stage_config: StageDefinition,
        manager: PipelineManager,
        current_params: Dict[str, Any],
    ):
        """
        Called by the main window to set the *context* for the viewer.
        This populates the dropdowns.
        """
        self.stage_config = stage_config
        self.manager = manager
        self.current_params = current_params

        # Enable the ribosome toggle now that we have a manager
        self.ribosome_toggle.setEnabled(True)

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

        self.info_text.setPlainText(
            f"Stage: {stage_config.name}\n"
            f"Found {len(current_stage_artifacts)} artifact files for this stage."
        )

        # 1. Add registered visualizations
        for viz_spec in stage_config.visualizations:
            self.viz_select_combo.addItem(viz_spec.name, userData=viz_spec)

        # 2. Add the "Raw File Inspector" as a fallback
        self.viz_select_combo.addItem(
            "Inspect Raw Artifact File", userData="RAW_FILE_INSPECTOR"
        )

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
            if ext == ".npy":
                return np.load(path, allow_pickle=True)
            elif ext in [".ply", ".stl", ".obj", ".vtk"]:
                return pv.read(path)
            elif ext == ".json":
                with open(path, "r") as f:
                    return json.load(f)
            elif ext == ".cif":
                return str(path)
            else:
                return str(path)
        except Exception as e:
            print(f"Error loading data from {path}: {e}")
            raise

    def execute_visualization(self, viz_spec: VisualizationSpec):
        """Gathers data and calls the registered viz function."""

        # Clear plotter but preserve ribosome overlay state
        ribosome_was_shown = self.show_ribosome_structure
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
                    kwargs[kwarg_name] = None  # Pass None to the function
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

            # Re-add ribosome overlay if it was shown
            if ribosome_was_shown:
                self.add_ribosome_overlay()

        except Exception as e:
            error_msg = f"Failed to execute visualization: {viz_spec.name}\n\n"
            error_msg += f"Error: {str(e)}\n"
            error_msg += f"{traceback.format_exc()}"
            self.info_text.setPlainText(error_msg)
            self.plotter.add_text(f"Error: {e}", position="upper_left", color="red")

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

        # Preserve ribosome overlay state
        ribosome_was_shown = self.show_ribosome_structure
        self.plotter.clear()
        self.plotter.enable_trackball_style()

        try:
            ext = artifact_path.suffix.lower()

            if ext in [".ply", ".stl", ".obj", ".vtk"]:
                mesh = pv.read(artifact_path)
                self.plotter.add_mesh(
                    mesh, color="lightblue", show_edges=True, opacity=0.8
                )
                info = f"Mesh: {artifact_path.name}\n"
                info += f"Points: {mesh.n_points}\n"
                info += f"Faces: {mesh.n_faces}\n"
                self.info_text.setPlainText(info)

            elif ext == ".npy":
                points = np.load(artifact_path)
                if len(points.shape) == 2 and points.shape[1] == 3:
                    cloud = pv.PolyData(points)
                    self.plotter.add_mesh(
                        cloud, color="red", point_size=5, render_points_as_spheres=True
                    )
                    info = f"Point Cloud: {artifact_path.name}\n"
                    info += f"Points: {len(points)}\n"
                    self.info_text.setPlainText(info)
                else:
                    self.info_text.setPlainText(
                        f"Numpy array: {artifact_path.name}\nShape: {points.shape}"
                    )

            elif ext == ".cif":
                info = f"Structure File: {artifact_path.name}\n"
                info += (
                    f"Size: {artifact_path.stat().st_size / (1024 * 1024):.2f} MB\n\n"
                )
                info += "Raw CIF preview not supported."
                self.info_text.setPlainText(info)
                self.plotter.add_text(
                    "CIF structure file (see info panel)", position="upper_left"
                )

            elif ext == ".json":
                with open(artifact_path) as f:
                    data = json.load(f)
                self.info_text.setPlainText(json.dumps(data, indent=2))
                self.plotter.add_text(
                    "JSON data (see info panel)", position="upper_left"
                )

            else:
                self.info_text.setPlainText(f"Cannot preview file type: {ext}")
                self.plotter.add_text(f"Cannot preview {ext}", position="upper_left")

            self.plotter.reset_camera()

            # Re-add ribosome overlay if it was shown
            if ribosome_was_shown:
                self.add_ribosome_overlay()

        except Exception as e:
            self.info_text.setPlainText(f"Error loading artifact:\n{str(e)}")
            print(traceback.format_exc())

    def export_artifact(self):
        """Export the currently selected raw artifact file."""
        current_index = self.raw_file_combo.currentIndex()
        if current_index < 0:
            return

        artifact_path: Path = self.raw_file_combo.itemData(current_index)
        if not artifact_path:
            return

        save_path, _ = QFileDialog.getSaveFileName(
            self,
            "Export Artifact",
            str(Path.home() / artifact_path.name),
            f"{artifact_path.suffix.upper()} files (*{artifact_path.suffix});;All files (*.*)",
        )

        if save_path:
            try:
                import shutil

                shutil.copy(artifact_path, save_path)
                self.info_text.append(f"\nExported to {save_path}")
            except Exception as e:
                self.info_text.append(f"\nExport failed: {e}")


def main():
    os.environ["QT_API"] = "pyqt5"

    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)

    app = QApplication(sys.argv)
    app.setStyle("Fusion")

    orchestrator = NPETPipelineOrchestrator()
    orchestrator.show()

    sys.exit(app.exec_())


if __name__ == "__main__":
    if "RIBETL_DATA" not in os.environ:
        print("ERROR: 'RIBETL_DATA' environment variable not set.")
        # You should probably exit here, but we'll let it try

    pv.set_plot_theme("paraview")
    pv.global_theme.multi_rendering_splitting_position = 0.5

    main()
