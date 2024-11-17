from dataclasses import dataclass
from typing import Dict, Set, Optional, Union, Type
from pathlib import Path
from loguru import logger
from pydantic import BaseModel

from .asset_types import AssetType, ModelT,ModelGenerator, RawGenerator, AssetDefinition

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
              AssetType.MMCIF                : "{asset_dir}/{pdb_id}.cif",
              AssetType.STRUCTURE_PROFILE    : "{asset_dir}/{pdb_id}.json",
              AssetType.PTC                  : "{asset_dir}/{pdb_id}_PTC.json",
            #   AssetType.THUMBNAIL            : "{asset_dir}/{pdb_id}.png",
            #   AssetType.CLASSIFICATION_REPORT: "{asset_dir}/classification_report_{pdb_id}.json",
            #   AssetType.NPET_MESH            : "{asset_dir}/TUNNELS/{pdb_id}_NPET_MESH.ply",
            # AssetType.RNA_HELICES          : "{asset_dir}/RNA_HELICES.json",
            # AssetType.TRNA_SITES           : "{asset_dir}/TRNA_SITES.json",
        }

        template = templates.get(asset_type)
        if not template:
            raise ValueError(f"No path template defined for asset type {asset_type}")

        return Path(template.format(asset_dir=asset_dir, pdb_id=pdb_id.upper()))

    def load_model(self, pdb_id: str, asset_type: AssetType) -> Optional[BaseModel]:
        """Load and validate a model from disk"""
        if not asset_type.requires_model():
            return None
            
        path = self.get_asset_path(pdb_id, asset_type)
        if not path.exists():
            raise FileNotFoundError(f"No {asset_type.name} asset found for {pdb_id}")
            
        model_cls = asset_type.model_type
        return model_cls.model_validate_json(path.read_text())

class RibosomeAssetManager:
    """Manages assets for ribosome structures"""
    def __init__(self, base_dir: Path):
        self.path_manager = AssetPathManager(base_dir)
        self._init_asset_definitions()

    def _init_asset_definitions(self):
        """Initialize asset definitions - now uses model information from AssetType"""
        self.assets: Dict[AssetType, AssetDefinition] = {}
        
        for asset_type in AssetType:
            self.assets[asset_type] = AssetDefinition(
                asset_type    = asset_type,
                dependencies  = asset_type.dependencies,
                path_template = "{pdb_id}",              # Basic template, actual paths handled by PathManager
                required      = True,                    # Could be made configurable per asset if needed
                generator     = None  # Will be set by registry
            )

    def register_generator(
        self, 
        asset_type: AssetType, 
        generator: Union[ModelGenerator[ModelT], RawGenerator]
    ) -> None:
        """Register a generator function for an asset type"""
        if asset_type not in self.assets:
            raise ValueError(f"Unknown asset type: {asset_type}")
        self.assets[asset_type].generator = generator

    async def verify_asset(self, pdb_id: str, asset_type: AssetType) -> bool:
        """Verify that an asset exists and optionally validate its model"""
        path = self.path_manager.get_asset_path(pdb_id, asset_type)
        if not path.exists():
            return False

        # For model-based assets, try to load and validate the model
        if asset_type.requires_model():
            try:
                self.path_manager.load_model(pdb_id, asset_type)
            except Exception as e:
                logger.warning(f"Asset exists but failed validation: {str(e)}")
                return False

        return True

    async def verify_all_assets(self, pdb_id: str) -> Dict[AssetType, bool]:
        """Verify all assets for a structure"""
        return {
            asset_type: await self.verify_asset(pdb_id, asset_type)
            for asset_type in AssetType
        }

    def load_dependency(self, pdb_id: str, asset_type: AssetType) -> Optional[BaseModel]:
        """Load a dependency model, returns None for raw assets"""
        return self.path_manager.load_model(pdb_id, asset_type)
