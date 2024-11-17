from dataclasses import dataclass
from enum import Enum
from typing import Protocol, Set, TypeVar, Optional, Any, Union
from pydantic import BaseModel
from pathlib import Path
from typing import TypeVar, Generic, Optional, Type, Tuple
from pydantic import BaseModel

from ribctl.lib.schema.types_ribosome import PTCInfo, RibosomeStructure

# Type for any pydantic model
ModelT = TypeVar('ModelT', bound=BaseModel)

class AssetInfo:
    def __init__(
        self, 
        name: str, 
        model: Optional[Type[BaseModel]] = None,
        dependencies: set[str] = set(),
        is_raw: bool = False  # New flag to explicitly mark raw assets
    ):
        self.name = name
        self.model = model
        self.dependencies = dependencies
        self.is_raw = is_raw

class AssetType(Enum):
    # Raw assets explicitly marked
    MMCIF     = AssetInfo("mmcif", model=None, is_raw=True)
    THUMBNAIL = AssetInfo("thumbnail", model=None, is_raw=True, dependencies={"STRUCTURE_PROFILE","MMCIF"})
    NPET_MESH = AssetInfo("npet_mesh", model=None, is_raw=True, dependencies={"STRUCTURE_PROFILE","MMCIF"})
    
    # Model-based assets
    STRUCTURE_PROFILE = AssetInfo("profile", model=RibosomeStructure)
    PTC               = AssetInfo("ptc", model=PTCInfo, dependencies={"STRUCTURE_PROFILE", "MMCIF"});
    
    @property
    def model_type(self) -> Optional[Type[BaseModel]]:
        return self.value.model

    @property
    def is_raw_asset(self) -> bool:
        return self.value.is_raw

    @property
    def dependencies(self) -> set['AssetType']:
        return {AssetType[d] for d in self.value.dependencies}

    def requires_model(self) -> bool:
        return self.model_type is not None


ModelT = TypeVar('ModelT', bound=BaseModel)

class ModelGenerator(Protocol[ModelT]):
    """Protocol for generators that produce model-based assets"""
    async def __call__(self, rcsb_id: str, dependencies: Optional[dict[str, Any]] = None) -> ModelT:
        ...

class RawGenerator(Protocol):
    """Protocol for generators that handle their own file I/O"""
    async def __call__(self, rcsb_id: str, output_path: Path) -> None:
        ...

@dataclass
class AssetDefinition(Generic[ModelT]):
    """Enhanced asset definition that knows about its model type"""
    asset_type    : AssetType
    dependencies  : Set[AssetType]
    path_template : str
    required      : bool
    generator     : Optional[Union[ModelGenerator[ModelT], RawGenerator]] = None
    schema_version: str = "0.1"

    @property
    def model_type(self) -> Optional[Type[BaseModel]]:
        return self.asset_type.model_type