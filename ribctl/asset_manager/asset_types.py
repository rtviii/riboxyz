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

class AssetInfo(Generic[ModelT]):
    """Associates an asset with its model type and other metadata"""
    def __init__(
        self, 
        name: str, 
        model: Optional[Type[ModelT]] = None,
        dependencies: set[str] = set()
    ):
        self.name = name
        self.model = model
        self.dependencies = dependencies

class AssetType(Enum):
    # Core Data Assets - note how each specifies its model type
    MMCIF                 = AssetInfo("mmcif", None)  # Raw file, no model
    STRUCTURE_PROFILE     = AssetInfo("profile", RibosomeStructure)
    # CLASSIFICATION_REPORT = AssetInfo("classification", ClassificationReport)

    # # Structural Analysis Assets
    PTC         = AssetInfo("ptc", PTCInfo, dependencies={"MMCIF", "STRUCTURE_PROFILE"})
    # NPET_MESH   = AssetInfo("npet", NPETMesh, dependencies={"MMCIF"})
    # RNA_HELICES = AssetInfo("rna_helices", RNAHelices)
    # TRNA_SITES  = AssetInfo("trna", TRNASites)

    # Visualization Assets
    THUMBNAIL = AssetInfo("thumbnail", None)  # Another raw file

    @property
    def model_type(self) -> Optional[Type[BaseModel]]:
        return self.value.model

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