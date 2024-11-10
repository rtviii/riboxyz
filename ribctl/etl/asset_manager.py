import asyncio
from dataclasses import dataclass
from enum import Enum, auto
from pprint import pprint
from typing import Optional, Set, Dict
from pathlib import Path
import json
from datetime import datetime
from ribctl import RIBETL_DATA


class AssetType(Enum):
    # Core Data Assets
    MMCIF = auto()
    STRUCTURE_PROFILE = auto()  # Structure profile
    CLASSIFICATION_REPORT = auto()

    # Structural Analysis Assets
    PTC = auto()  # PTC site analysis
    NPET_MESH = auto()  # Tunnel analysis
    RNA_HELICES = auto()  # RNA helix annotations
    TRNA_SITES = auto()  # tRNA binding sites

    # Visualization Assets
    THUMBNAIL = auto()  # Structure thumbnail


@dataclass
class AssetDefinition:
    """Defines an asset and its properties"""

    asset_type: AssetType
    dependencies: Set[AssetType]
    path_template: str
    required: bool
    schema_version: str = "0.1"


class AssetPathManager:
    """Manages asset paths and directory structure"""

    def __init__(self, base_dir: Path):
        self.base_dir = Path(base_dir)

    def get_asset_dir(self, pdb_id: str) -> Path:
        """Get the base directory for a structure's assets"""
        return self.base_dir / pdb_id.upper()

    def get_asset_path(self, pdb_id: str, asset_type: AssetType) -> Path:
        """Get the path for a specific asset"""
        asset_dir = self.get_asset_dir(pdb_id)

        # Define path templates for each asset type
        templates = {
            AssetType.MMCIF: "{asset_dir}/{pdb_id}.cif",
            AssetType.STRUCTURE_PROFILE: "{asset_dir}/{pdb_id}.json",
            AssetType.PTC: "{asset_dir}/{pdb_id}_PTC.json",
            AssetType.THUMBNAIL: "{asset_dir}/{pdb_id}.png",
            AssetType.CLASSIFICATION_REPORT: "{asset_dir}/classification_report_{pdb_id}.json",
            AssetType.NPET_MESH: "{asset_dir}/TUNNELS",
            AssetType.RNA_HELICES: "{asset_dir}/RNA_HELICES.json",
            AssetType.TRNA_SITES: "{asset_dir}/TRNA_SITES.json",
        }

        template = templates.get(asset_type)
        if not template:
            raise ValueError(f"No path template defined for asset type {asset_type}")

        return Path(template.format(asset_dir=asset_dir, pdb_id=pdb_id.upper()))


class RibosomeAssetManager:
    """Manages assets for ribosome structures"""

    def __init__(self, base_dir: Path):
        self.path_manager = AssetPathManager(base_dir)
        self._init_asset_definitions()

    def _init_asset_definitions(self):
        """Initialize asset definitions with their dependencies"""
        self.assets: Dict[AssetType, AssetDefinition] = {
            AssetType.MMCIF: AssetDefinition(
                asset_type=AssetType.MMCIF,
                dependencies=set(),  # No dependencies
                path_template="{pdb_id}.cif",
                required=True,
            ),
            AssetType.STRUCTURE_PROFILE: AssetDefinition(
                asset_type=AssetType.STRUCTURE_PROFILE,
                dependencies=set(),
                path_template="{pdb_id}.json",
                required=True,
            ),
            AssetType.PTC: AssetDefinition(
                asset_type=AssetType.PTC,
                dependencies={AssetType.MMCIF, AssetType.STRUCTURE_PROFILE},
                path_template="{pdb_id}_PTC.json",
                required=True,
            ),
        }

    async def verify_asset(self, pdb_id: str, asset_type: AssetType) -> bool:
        """Verify that an asset exists and is valid"""
        path = self.path_manager.get_asset_path(pdb_id, asset_type)
        if not path.exists():
            return False

        return True

    async def verify_all_assets(self, pdb_id: str) -> Dict[AssetType, bool]:
        """Verify all assets for a structure"""
        return { asset_type: await self.verify_asset(pdb_id, asset_type) for asset_type in AssetType }


async def example_usage():

    manager      = RibosomeAssetManager(RIBETL_DATA)
    pdb_id       = "4UG0"
    asset_status = await manager.verify_all_assets(pdb_id)
    pprint(asset_status)

    # Get path for specific asset
    ptc_path = manager.path_manager.get_asset_path(pdb_id, AssetType.PTC)

    # Check dependencies for PTC asset
    ptc_definition = manager.assets[AssetType.PTC]
    print(f"PTC asset dependencies: {ptc_definition.dependencies}")


