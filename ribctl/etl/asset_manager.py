from dataclasses import dataclass
from enum import Enum, auto
from typing import Dict, Set, Optional, Protocol
from pathlib import Path
from loguru import logger
from pathlib import Path
import functools
from loguru import logger
from typing import Dict, Optional


class AssetGenerator(Protocol):
    """Protocol that all asset generators must implement"""

    async def __call__(self, rcsb_id: str, overwrite: bool = False) -> None: ...


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
    generator: Optional[AssetGenerator] = None  # Only new field added
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

    def register_generator(
        self, asset_type: AssetType, generator: AssetGenerator
    ) -> None:
        """One new method to register generators"""
        if asset_type not in self.assets:
            raise ValueError(f"Unknown asset type: {asset_type}")
        self.assets[asset_type].generator = generator

    async def verify_asset(self, pdb_id: str, asset_type: AssetType) -> bool:
        """Verify that an asset exists and is valid"""
        path = self.path_manager.get_asset_path(pdb_id, asset_type)
        if not path.exists():
            return False
        return True

    async def verify_all_assets(self, pdb_id: str) -> Dict[AssetType, bool]:
        """Verify all assets for a structure"""
        return {
            asset_type: await self.verify_asset(pdb_id, asset_type)
            for asset_type in AssetType
        }


class AssetGenRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager
        self.path_manager = manager.path_manager

    def register(self, asset_type: AssetType):
        """Decorator to register a generator function"""

        def decorator(func: AssetGenerator) -> AssetGenerator:
            @functools.wraps(func)
            async def wrapped(rcsb_id: str, overwrite: bool = False) -> None:
                output_path = self.path_manager.get_asset_path(rcsb_id, asset_type)
                try:
                    # Check if asset exists
                    if output_path.exists():
                        if not overwrite:
                            logger.info(
                                f"Asset exists at {output_path}, overwrite=False, skipping"
                            )
                            return
                        logger.info(
                            f"Asset exists at {output_path}, overwrite=True, regenerating"
                        )

                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    await func(rcsb_id, overwrite)
                    logger.success(f"Generated {func.__name__} for {rcsb_id}")

                except Exception as e:
                    logger.exception(f"Failed {func.__name__} for {rcsb_id}: {str(e)}")
                    raise

            # Register the wrapped function
            self.manager.register_generator(asset_type, wrapped)
            return wrapped

        return decorator

    async def generate_asset(
        self, rcsb_id: str, asset_type: AssetType, force: bool = False
    ) -> None:
        """Generate a single asset and its dependencies"""
        logger.info(f"Generating {asset_type.name} for {rcsb_id}")

        asset_def = self.manager.assets[asset_type]

        # Handle dependencies first
        if asset_def.dependencies:
            logger.info(f"Processing dependencies for {asset_type.name}")
            for dep in asset_def.dependencies:
                await self.generate_asset(rcsb_id, dep, force)

        # Generate the asset
        if not asset_def.generator:
            raise ValueError(f"No generator registered for {asset_type}")

        await asset_def.generator(rcsb_id, force)

    async def generate_multiple(
        self, rcsb_id: str, asset_types: list[AssetType], force: bool = False
    ) -> None:
        """Generate multiple assets for a structure"""
        for asset_type in asset_types:
            await self.generate_asset(rcsb_id, asset_type, force)
