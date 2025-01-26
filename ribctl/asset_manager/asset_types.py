from dataclasses import dataclass
from enum import Enum
from typing import Protocol, Set, TypeVar, Optional, Any, Union
from pydantic import BaseModel
from pathlib import Path
from typing import TypeVar, Generic, Type
from ribctl.lib.schema.types_ribosome import (
    ConstrictionSite,
    PTCInfo,
    RibosomeStructure,
)
from ribctl import RIBETL_DATA

# Type for any pydantic model
ModelT = TypeVar("ModelT", bound=BaseModel)


class AssetInfo:
    def __init__(
        self,
        name: str,
        path_template: str,
        model: Optional[Type[BaseModel]] = None,
        dependencies: set[str] = set(),
        is_raw: bool = False,
    ):
        self.name = name
        self.path_template = path_template
        self.model = model
        self.dependencies = dependencies
        self.is_raw = is_raw

class PathResolver:
    _base_dir = Path(RIBETL_DATA)

    @classmethod
    def set_base_dir(cls, path: Path):
        cls._base_dir = Path(path)

    @classmethod
    def get_base_dir(cls) -> Path:
        return cls._base_dir

class AssetType(Enum):

    @classmethod
    def get_base_dir(cls) -> Path:
        return Path(cls._base_dir)

    def get_path(self, pdb_id: str) -> Path:
        """Get path for this asset type"""
        asset_dir = PathResolver.get_base_dir() / pdb_id.upper()
        return Path(self.value.path_template.format(
            asset_dir=asset_dir,
            pdb_id=pdb_id.upper()
        ))

    # Raw assets explicitly marked
    MMCIF = AssetInfo(
        "mmcif", path_template="{asset_dir}/{pdb_id}.cif", model=None, is_raw=True
    )

    STRUCTURE_PROFILE = AssetInfo(
        "profile", path_template="{asset_dir}/{pdb_id}.json", model=RibosomeStructure
    )

    CONSTRICTION_SITE = AssetInfo(
        "constriction_site",
        path_template="{asset_dir}/{pdb_id}_CONSTRICTION_SITE.json",
        model=ConstrictionSite,
        dependencies={"STRUCTURE_PROFILE", "MMCIF"},
    )

    PTC = AssetInfo(
        "ptc",
        path_template="{asset_dir}/{pdb_id}_PTC.json",
        model=PTCInfo,
        dependencies={"STRUCTURE_PROFILE", "MMCIF"},
    )

    THUMBNAIL = AssetInfo(
        "thumbnail",
        path_template="{asset_dir}/{pdb_id}.png",
        model=None,
        is_raw=True,
        dependencies={"STRUCTURE_PROFILE", "MMCIF"},
    )

    NPET_MESH = AssetInfo(
        "npet_mesh",
        path_template="{asset_dir}/TUNNELS/{pdb_id}_NPET_MESH.ply",
        model=None,
        is_raw=True,
        dependencies={
            "STRUCTURE_PROFILE",
            "MMCIF",
            "ALPHA_SHAPE",
            "CONSTRICTION_SITE",
            "PTC",
        },
    )

    ALPHA_SHAPE = AssetInfo(
        "alpha_shape",
        path_template="{asset_dir}/{pdb_id}_ALPHA_SHAPE.ply",
        model=None,
        dependencies={"MMCIF"},
    )

    SRL = AssetInfo(
        "sarcin_ricin_loop",
        path_template="{asset_dir}/{pdb_id}_SRL.json",
        model=None,
        is_raw=True,
        dependencies={"MMCIF", "STRUCTURE_PROFILE"},
    )

    L7L12STALK = AssetInfo(
        "l7l12_stalk",
        path_template="{asset_dir}/{pdb_id}_L7L12STALK.json",
        model=None,
        is_raw=True,
        dependencies={"MMCIF", "STRUCTURE_PROFILE"},
    )

    PSTALK = AssetInfo(
        "p_stalk",
        path_template="{asset_dir}/{pdb_id}_PSTALK.json",
        model=None,
        is_raw=True,
        dependencies={"MMCIF", "STRUCTURE_PROFILE"},
    )

    @property
    def model_type(self) -> Optional[Type[BaseModel]]:
        return self.value.model

    @property
    def is_raw_asset(self) -> bool:
        return self.value.is_raw

    @property
    def dependencies(self) -> set["AssetType"]:
        return {AssetType[d] for d in self.value.dependencies}

    def requires_model(self) -> bool:
        return self.model_type is not None


class ModelGenerator(Protocol[ModelT]):
    """Protocol for generators that produce model-based assets"""

    async def __call__(
        self, rcsb_id: str, dependencies: Optional[dict[str, Any]] = None
    ) -> ModelT: ...


class RawGenerator(Protocol):
    """Protocol for generators that handle their own file I/O"""

    async def __call__(self, rcsb_id: str, output_path: Path) -> None: ...


@dataclass
class AssetDefinition(Generic[ModelT]):
    """Enhanced asset definition that knows about its model type"""

    asset_type: AssetType
    dependencies: Set[AssetType]
    required: bool
    generator: Optional[Union[ModelGenerator[ModelT], RawGenerator]] = None
    schema_version: str = "0.1"

    @property
    def model_type(self) -> Optional[Type[BaseModel]]:
        return self.asset_type.model_type
