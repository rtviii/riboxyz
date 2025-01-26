# Asset Management System Documentation

## Overview
This system manages various types of assets (files, data structures) associated with molecular structures, handling both raw files and typed (Pydantic model) assets.

## Key Components

### AssetType Enum
- Defines available asset types using `AssetInfo` class
- Specifies whether asset is raw or typed (Pydantic model)
- Declares dependencies between assets
- Example: `MMCIF = AssetInfo("mmcif", model=None, is_raw=True)`

### Asset Registry
- Manages typed assets using Pydantic models
- Provides decorator for registering asset generators
- Handles dependency resolution
```python
@main_registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    return await ETLCollector(rcsb_id).generate_profile()
```

### Raw Asset Handler
- Manages file-based assets without Pydantic models
- Handlers registered via `register_handler()`
- Built-in handlers defined in `_register_default_handlers()`

### Parallel Acquisition
- Processes multiple structures concurrently
- Respects asset dependencies
- Controls concurrency via semaphores
- Tracks success/failure via `AcquisitionResult`

## CLI Commands
- `get`: Process specific structure(s)
- `get-all`: Process all structures in parallel chunks

## Adding New Assets
1. Add entry to `AssetType` enum
2. For typed assets: Register generator with `@main_registry.register`
3. For raw assets: Add handler to `RawAssetHandler`