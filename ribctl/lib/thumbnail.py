import os
import subprocess
from pathlib import Path
from typing import Optional
import logging
from dataclasses import dataclass

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Environment variables
CHIMERAX_BINARY_PATH = os.getenv("CHIMERAX_BINARY_PATH")
RIBXZ_ROOT = os.getenv("RIBXZ_ROOT")
RIBETL_DATA = os.getenv("RIBETL_DATA")


@dataclass
class WorkflowConfig:
    """Configuration for the visualization workflow"""
    chimerax_path  : str
    project_root   : Path
    output_dir     : Path
    fps            : int = 30
    gif_width      : int = 720
    gif_crop_factor: float = 0.5
    gif_colors     : int = 256


class RibosomeVisualizer:
    def __init__(self, config: WorkflowConfig):
        self.config = config
        self._validate_config()
        
    def _validate_config(self):
        """Validate the configuration and create directories if needed"""
        if not self.config.chimerax_path:
            self.config.chimerax_path = CHIMERAX_BINARY_PATH
            print(self.config.chimerax_path)
            if not self.config.chimerax_path:
                raise ValueError("CHIMERAX_PATH environment variable not set")
                
        # Create output directory if it doesn't exist
        os.makedirs(self.config.output_dir, exist_ok=True)
        
    def _get_script_path(self, script_name: str) -> str:
        """Get the absolute path for a ChimeraX script relative to project root"""
        script_path = os.path.join(self.config.project_root,'ribctl', "lib", "chimerax", script_name)
        if not os.path.exists(script_path):
            raise FileNotFoundError(f"Script not found: {script_path}")
        return script_path

    def produce_visualization(self, rcsb_id: str, output_type: str = "all") -> dict:
        """
        Produce visualization for given RCSB ID
        output_type: 'mp4', 'gif', 'png', or 'all'
        Returns dict with paths to generated files
        """
        output_files = {}
        
        # Generate MP4
        mp4_path = os.path.join(self.config.output_dir, f"{rcsb_id}.mp4")
        self._render_movie(rcsb_id, mp4_path)
        output_files["mp4"] = mp4_path
        
        if output_type in ["gif", "all"]:
            gif_path = os.path.join(self.config.output_dir, f"{rcsb_id}.gif")
            self._convert_to_gif(mp4_path, gif_path)
            output_files["gif"] = gif_path
            
        if output_type in ["png", "all"]:
            png_path = os.path.join(self.config.output_dir, f"{rcsb_id}.png")
            self._extract_frame(mp4_path, png_path)
            output_files["png"] = png_path
            
        return output_files

    def _render_movie(self, rcsb_id: str, output_path: str):
        """Render movie using ChimeraX"""
        cmd_scripts = [
            self._get_script_path("cmd_ribetl.py"),
            self._get_script_path("cmd_ribrepr.py"),
            self._get_script_path("cmd_ribmovie.py")
        ]
        
        # Split the ChimeraX command into components (e.g., "flatpak run edu.ucsf.rbvi.ChimeraX" -> ["flatpak", "run", "edu.ucsf.rbvi.ChimeraX"])
        chimerax_cmd = self.config.chimerax_path.split() + ["--nogui", "--offscreen", "--cmd"]
        
        # Build ChimeraX command
        script_loads = "; ".join(f"open {script}" for script in cmd_scripts)
        full_cmd = f"{script_loads}; ribmovie {rcsb_id}; save {output_path}; close all"
        
        try:
            subprocess.run([*chimerax_cmd, full_cmd], check=True)
            logger.info(f"Successfully rendered movie for {rcsb_id}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to render movie: {e}")
            raise

    def _convert_to_gif(self, input_path: str, output_path: str):
        """Convert MP4 to GIF with high quality settings"""
        filter_complex = (
            f"fps={self.config.fps},"
            f"scale={self.config.gif_width}:-1:flags=bicubic,"
            f"crop=iw*{self.config.gif_crop_factor}:ih:iw*{(1-self.config.gif_crop_factor)/2}:0,"
            f"split[s0][s1];"
            f"[s0]palettegen=max_colors={self.config.gif_colors}:reserve_transparent=0:stats_mode=full[p];"
            f"[s1][p]paletteuse=dither=none"
        )
        
        cmd = [
            "ffmpeg",
            "-y",
            "-i",
            str(input_path),
            "-filter_complex",
            filter_complex,
            "-loop",
            "0",
            str(output_path)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"Successfully converted {input_path} to GIF")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to convert to GIF: {e}")
            raise

    def _extract_frame(self, input_path: str, output_path: str, timestamp: float = 0):
        """Extract a single frame from the video"""
        cmd = [
            "ffmpeg",
            "-y",
            "-i",
            str(input_path),
            "-ss",
            str(timestamp),
            "-vframes",
            "1",
            "-vf",
            f"scale={self.config.gif_width}:-1:flags=bicubic,"
            f"crop=iw*{self.config.gif_crop_factor}:ih:iw*{(1-self.config.gif_crop_factor)/2}:0",
            "-compression_level",
            "0",
            str(output_path)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            logger.info(f"Successfully extracted frame to {output_path}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to extract frame: {e}")
            raise


# Example usage
if __name__ == "__main__":
    # Configuration can be loaded from environment variables or config file
    rcsb_id = "7K00"
    config = WorkflowConfig(
        chimerax_path=CHIMERAX_BINARY_PATH,
        project_root=RIBXZ_ROOT,
        output_dir=os.path.join(RIBETL_DATA, rcsb_id.upper())
    )
    
    visualizer = RibosomeVisualizer(config)
    try:
        output_files = visualizer.produce_visualization(rcsb_id)
        logger.info(f"Generated files: {output_files}")
    except Exception as e:
        logger.error(f"Failed to process {rcsb_id}: {e}")