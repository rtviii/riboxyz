Assuming you're in the riboxyz repo root:

## Riboxyz mode (you have local assets)

Single structure:
```bash
python -m ribctl.lib.npet2 run 7K00
```

Multiple structures:
```bash
python -m ribctl.lib.npet2 run 7K00 4UG0 3J7Z
```

Parallel:
```bash
python -m ribctl.lib.npet2 run 7K00 4UG0 3J7Z 5AFI -j 4
```

From a file (one ID per line, `#` comments allowed):
```bash
python -m ribctl.lib.npet2 run --from-file structures.txt -j 8
```

Custom output directory:
```bash
python -m ribctl.lib.npet2 run 7K00 --output-dir ./my_runs
```

## Config overrides

Tweak cylinder geometry:
```bash
python -m ribctl.lib.npet2 run 7K00 --cylinder-radius 40 --cylinder-height 130
```

Change voxel resolution:
```bash
python -m ribctl.lib.npet2 run 7K00 --voxel-size 0.5
```

Tune DBSCAN:
```bash
python -m ribctl.lib.npet2 run 7K00 --dbscan-coarse-eps 6.0 --dbscan-coarse-min-samples 400
```

Skip meshing (faster, just get the point clouds):
```bash
python -m ribctl.lib.npet2 run 7K00 --no-mesh
```

Full config from JSON:
```bash
python -m ribctl.lib.npet2 run 7K00 --config-json my_config.json
```

To see what the default config looks like:
```bash
python -m ribctl.lib.npet2 show-config
```

That dumps the full `RunConfig` as JSON, which you can save, edit, and pass back with `--config-json`.

## Standalone mode (no riboxyz repo)

This is for when someone has just an mmCIF file and gets profile/landmarks from the API or local JSON files:

```bash
# Fetch profile + landmarks from the riboxyz API
python -m ribctl.lib.npet2 run 7K00 \
    --mode standalone \
    --mmcif /data/7K00.cif \
    --api-url http://riboxyz.server:8000

# Everything from local files
python -m ribctl.lib.npet2 run 7K00 \
    --mode standalone \
    --mmcif /data/7K00.cif \
    --profile /data/7K00_profile.json \
    --landmarks /data/7K00_landmarks.json
```

For batch standalone with a directory of files:
```bash
# Expects 7K00.cif, 4UG0.cif, etc. in /data/mmcif/
# Expects 7K00_profile.json, etc. in /data/profiles/
# Expects 7K00_landmarks.json, etc. in /data/landmarks/
python -m ribctl.lib.npet2 run 7K00 4UG0 \
    --mode standalone \
    --mmcif /data/mmcif/ \
    --profile /data/profiles/ \
    --landmarks /data/landmarks/ \
    -j 4
```

The landmarks JSON format is:
```json
{
    "ptc": {"location": [x, y, z]},
    "constriction": {"location": [x, y, z]}
}
```

## Environment variables

Instead of CLI flags, you can set:
```bash
export NPET2_ROOT=~/npet2_data
export NPET2_RUNS_ROOT=~/npet2_data/runs
export NPET2_POISSON_RECON_BIN=/usr/local/bin/PoissonRecon
export NPET2_RIBOXYZ_API_URL=http://localhost:8000
```