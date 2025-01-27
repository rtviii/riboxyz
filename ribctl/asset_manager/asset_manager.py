from typing import Dict, Set, Optional, Union, Type
from pathlib import Path
from loguru import logger
from pydantic import BaseModel
from .asset_types import (
    AssetType,
    ModelT,
    ModelGenerator,
    PathResolver,
    RawGenerator,
    AssetDefinition,
)

class RibosomeAssetManager:
    """Manages assets for ribosome structures, handling generation, validation and loading"""

    def __init__(self, base_dir: Optional[Path] = None):
        """
        Initialize the asset manager
        Args:
            base_dir: Optional override for the base directory. If not provided, uses RIBETL_DATA
        """
        if base_dir:
            PathResolver.set_base_dir(base_dir)
        self._init_asset_definitions()

    def _init_asset_definitions(self):
        """Initialize asset definitions for all asset types"""
        self.assets: Dict[AssetType, AssetDefinition] = {}
        
        for asset_type in AssetType:
            self.assets[asset_type] = AssetDefinition(
                asset_type=asset_type,
                dependencies=asset_type.dependencies,
                required=True,  # Could be made configurable per asset if needed
                generator=None  # Will be set by registry
            )

    def register_generator(
        self,
        asset_type: AssetType,
        generator: Union[ModelGenerator[ModelT], RawGenerator],
    ) -> None:
        """
        Register a generator function for an asset type
        
        Args:
            asset_type: The asset type to register a generator for
            generator: The generator function to register
        Raises:
            ValueError: If the asset type is unknown
        """
        if asset_type not in self.assets:
            raise ValueError(f"Unknown asset type: {asset_type}")
        self.assets[asset_type].generator = generator

    def load_model(self, pdb_id: str, asset_type: AssetType) -> Optional[BaseModel]:
        """
        Load and validate a model from disk
        
        Args:
            pdb_id: The PDB ID of the structure
            asset_type: The type of asset to load
        Returns:
            The loaded model if successful, None if asset doesn't require a model
        Raises:
            FileNotFoundError: If the asset file doesn't exist
        """
        if not asset_type.requires_model():
            return None

        path = asset_type.get_path(pdb_id)
        if not path.exists():
            raise FileNotFoundError(f"No {asset_type.name} asset found for {pdb_id}")

        model_cls = asset_type.model_type
        return model_cls.model_validate_json(path.read_text())

    async def verify_asset(self, pdb_id: str, asset_type: AssetType) -> bool:
        """
        Verify that an asset exists and optionally validate its model
        
        Args:
            pdb_id: The PDB ID of the structure
            asset_type: The type of asset to verify
        Returns:
            True if asset exists and is valid, False otherwise
        """
        path = asset_type.get_path(pdb_id)
        if not path.exists():
            return False

        # For model-based assets, try to load and validate the model
        if asset_type.requires_model():
            try:
                self.load_model(pdb_id, asset_type)
            except Exception as e:
                logger.warning(f"Asset exists but failed validation: {str(e)}")
                return False

        return True

    async def verify_all_assets(self, pdb_id: str) -> Dict[AssetType, bool]:
        """
        Verify all assets for a structure
        
        Args:
            pdb_id: The PDB ID of the structure
        Returns:
            Dictionary mapping asset types to their verification status
        """
        return {
            asset_type: await self.verify_asset(pdb_id, asset_type)
            for asset_type in AssetType
        }

    def load_dependency(
        self, pdb_id: str, asset_type: AssetType
    ) -> Optional[BaseModel]:
        """
        Load a dependency model
        
        Args:
            pdb_id: The PDB ID of the structure
            asset_type: The type of asset to load
        Returns:
            The loaded model if it requires one, None for raw assets
        """
        return self.load_model(pdb_id, asset_type)
