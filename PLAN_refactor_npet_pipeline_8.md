Alrighty, i have an NPET extraction pipeline (that i have just ported out of my bigger riboxyz) codebase that i want to separate from outer repo in which it is embedded (riboxyz). The tricky bit is that it is somewhat reliant on the riboxyz code and interfaces for now, but only weakly so. My intention is to package npet2 code separately and share it as a docker container and without any reliance on the `riboxyz` code at all so let's get rid of all the dependencies between the two. I think most of the things in `kdtree` and `alphalib` can be directly copied/ported. The interfaces for ptc and constriction sites we can talk about.

Here is the full layout of the repo:

```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  tree -L 6 -I 'node_modules|venv|__pycache__|profiles|cache|debug_output|*.fasta|*.csv|assets_*|staticfiles|api|assets|*.png|TUBETL_DATA|*.pkl|*hmm|*fasta|npet|*.mdx|*.ts.map|*.d.ts|nightingale|NPET2' .  e
.
├── __scripts
│   ├── gettrna.py
│   ├── ligclass_distill_composite.py
│   ├── ligclass_distill.py
│   ├── ligclass.py
│   ├── maps_acquire.py
│   ├── masif_plugin.py
│   ├── merge_lig_info.py
│   ├── mesh_grid.py
│   ├── move_classification.sh
│   ├── move_tunnels.sh
│   ├── narstructs_tally
│   ├── npet2_slice_viewer.py
│   ├── npet2_view.py
│   ├── npet2_viz_enhanced.py
│   ├── npet2_viz_run.py
│   ├── npy_to_pdb.py
│   ├── process_prediction_7k00.py
│   ├── pymol_visualtion.py
│   ├── rcsb_cumulative_entries_barplot.py
│   ├── reclassify.py
│   ├── ribxz_chimerax
│   │   ├── build
│   │   │   ├── bdist.macosx-10.9-universal2
│   │   │   └── lib
│   │   │       └── ribxz_chimerax
│   │   │           ├── __init__.py
│   │   │           ├── cmd_loci.py
│   │   │           ├── cmd_polymers.py
│   │   │           ├── cmd_registry.py
│   │   │           └── io.py
│   │   ├── bundle_info.xml
│   │   ├── dist
│   │   │   └── ribxz_chimerax-0.1-py3-none-any.whl
│   │   ├── ribxz_chimerax.egg-info
│   │   │   ├── dependency_links.txt
│   │   │   ├── PKG-INFO
│   │   │   ├── requires.txt
│   │   │   ├── SOURCES.txt
│   │   │   └── top_level.txt
│   │   └── src
│   │       ├── __init__.py
│   │       ├── cmd_loci.py
│   │       ├── cmd_polymers.py
│   │       ├── cmd_registry.py
│   │       └── io.py
│   ├── stats.json
│   ├── uniprot_seeds_query_record.py
│   └── williamson_assembly.py
├── dirtocontext.py
├── docker-compose.yml
├── Dockerfile-django
├── docs.md
├── muscle3.8.31_src.tar.gz
├── ncbi_taxonomy.sqlite
├── neo4j_ribosome
│   ├── __archive
│   │   ├── cypher_ops
│   │   │   ├── cypher_exec
│   │   │   ├── neo4j_commit_structure.sh
│   │   │   └── neo4j_seed_db_ontology.sh
│   │   └── riboxyz_seed_data
│   │       ├── ban-pfam-map-lsu.json
│   │       ├── ban-pfam-map-ssu.json
│   │       ├── interpro-base.json
│   │       ├── interpro-go-1.json
│   │       ├── interpro-go-2.json
│   │       ├── interpro-go-3.json
│   │       ├── interpro-go-4.json
│   │       ├── package-lock.json
│   │       ├── package.json
│   │       ├── pfam-interpro-1.json
│   │       ├── pfam-interpro-2.json
│   │       ├── pfam-interpro-3.json
│   │       ├── pfam-interpro-4.json
│   │       └── pfam-to-interpro-map.json
│   ├── __init__.py
│   ├── cypher
│   │   └── list_filtered_structs.cypher
│   ├── db_driver.py
│   ├── db_lib_builder.py
│   ├── db_lib_reader.py
│   ├── node_ligand.py
│   ├── node_phylogeny.py
│   ├── node_polymer.py
│   ├── node_protein.py
│   ├── node_rna.py
│   └── node_structure.py
├── neo4j.conf
├── neo4j.old.conf
├── notes
│   ├── binding_affinities.md
│   ├── docs
│   │   ├── architecture_environment_configs.md
│   │   └── update.md
│   ├── factors.md
│   ├── functional_sites.md
│   ├── general_directions.md
│   ├── hmm-based-classification.md
│   └── pymol.md
├── NPET_cli.md
├── npet_orchestrator.py
├── NPET_README.md
├── npet2_viewer_usage_examples.md
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_6.md
├── PLAN_refactor_npet_pipeline_7_cleanup.md
├── q_mass_runs_and_packaging.md
├── q_poisson_still_fails.md
├── q.md
├── ribctl
│   ├── __init__.py
│   ├── asset_manager
│   │   ├── asset_manager.py
│   │   ├── asset_registry.py
│   │   ├── asset_types.py
│   │   ├── doc.md
│   │   └── parallel_acquisition.py
│   ├── classifyre.json
│   ├── etl
│   │   ├── __init__.py
│   │   ├── etl_collector.py
│   │   └── gql_querystrings.py
│   ├── global_ops.py
│   ├── lib
│   │   ├── __libseq.py
│   │   ├── chimerax
│   │   │   ├── _cmd_ribrepr.py
│   │   │   ├── cmd_chainsplitter.py
│   │   │   ├── cmd_ligvis.py
│   │   │   ├── cmd_ribetl.py
│   │   │   ├── cmd_ribmovie.py
│   │   │   ├── cmd_ribrepr.py
│   │   │   ├── cmds_all.py
│   │   │   ├── ffmpeg_convert.sh
│   │   │   ├── ffmpeg_firstframe.sh
│   │   │   ├── gen_movies.py
│   │   │   ├── loop_ligvis.py
│   │   │   ├── loop_movies.py
│   │   │   ├── loop_split_chains.py
│   │   │   ├── notes.md
│   │   │   ├── produce_gif.py
│   │   │   └── thumbnails_from_ribetl.sh
│   │   ├── enumunion.py
│   │   ├── info.py
│   │   ├── landmarks
│   │   │   ├── __ptc_via_doris.py
│   │   │   ├── constriction_site.py
│   │   │   ├── notes.md
│   │   │   ├── ptc_via_trna.py
│   │   │   └── rrna_helices
│   │   │       ├── convert.py
│   │   │       ├── ecoli_7K00.json
│   │   │       ├── rrna_helices.py
│   │   │       ├── thermus_1VY4.json
│   │   │       └── yeast_4V88.json
│   │   ├── libbsite.py
│   │   ├── libhmm.py
│   │   ├── libmsa.py
│   │   ├── libseq.py
│   │   ├── libtax.py
│   │   ├── npet2
│   │   │   ├── __init__.py
│   │   │   ├── __main__.py
│   │   │   ├── adapters
│   │   │   │   ├── riboxyz_providers.py
│   │   │   │   └── standalone_providers.py
│   │   │   ├── backends
│   │   │   │   ├── __init__.py
│   │   │   │   ├── clustering_io.py
│   │   │   │   ├── grid_occupancy.py
│   │   │   │   ├── legacy
│   │   │   │   └── meshing.py
│   │   │   ├── core
│   │   │   │   ├── cache.py
│   │   │   │   ├── config.py
│   │   │   │   ├── interfaces.py
│   │   │   │   ├── manifest.py
│   │   │   │   ├── pipeline.py
│   │   │   │   ├── polymer_enum.py
│   │   │   │   ├── ribosome_types.py
│   │   │   │   ├── run_id.py
│   │   │   │   ├── settings.py
│   │   │   │   ├── store.py
│   │   │   │   ├── structure_selection.py
│   │   │   │   └── types.py
│   │   │   ├── run.py
│   │   │   └── stages
│   │   │       ├── bootstrap.py
│   │   │       ├── grid_refine.py
│   │   │       └── legacy_minimal.py
│   │   ├── nsearch_gemmi.py
│   │   ├── ribosome_types
│   │   ├── schema
│   │   │   ├── __init__.py
│   │   │   ├── primitives.py
│   │   │   ├── types_binding_site.py
│   │   │   └── types_ribosome.py
│   │   ├── seq_project_many_to_one.py
│   │   ├── thumbnail.py
│   │   ├── types
│   │   │   └── polymer
│   │   │       ├── __init__.py
│   │   │       ├── base.py
│   │   │       ├── hierarchies.py
│   │   │       └── types.py
│   │   └── utils.py
│   ├── logger_config.py
│   ├── logs
│   │   ├── etl.log
│   │   └── loggers.py
│   ├── ribd.py
│   └── ribosome_ops.py
├── taxdump.tar.gz
└── test_npet2.py
```

And now let me show you the actual code codebase.


ribctl/lib/npet2/adapters/riboxyz_providers.py
```py
# ribctl/lib/npet2/adapters/riboxyz_providers.py
"""
Providers that bridge riboxyz internals -> npet2 interfaces.

These are only usable inside the riboxyz repo where RibosomeOps, AssetType,
PTC_location, get_constriction are available. For standalone use, see
standalone_providers.py.
"""
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from ribctl.lib.npet2.core.ribosome_types import (
    RibosomeProfile,
    ProteinEntry,
    RNAEntry,
    PolymerEntry,
    AssemblyInfo,
    AssemblyPolymerInstance,
    AssemblyNonpolymerInstance,
)


def _convert_profile(ribo_struct) -> RibosomeProfile:
    """Convert a riboxyz RibosomeStructure to the npet2-internal RibosomeProfile."""

    def _nomenclature_strings(poly) -> list[str]:
        return [n.value if hasattr(n, "value") else str(n) for n in (poly.nomenclature or [])]

    proteins = [
        ProteinEntry(
            auth_asym_id=p.auth_asym_id,
            assembly_id=p.assembly_id,
            nomenclature=_nomenclature_strings(p),
            rcsb_pdbx_description=p.rcsb_pdbx_description,
            entity_poly_seq_length=p.entity_poly_seq_length,
        )
        for p in (ribo_struct.proteins or [])
    ]
    rnas = [
        RNAEntry(
            auth_asym_id=r.auth_asym_id,
            assembly_id=r.assembly_id,
            nomenclature=_nomenclature_strings(r),
            rcsb_pdbx_description=r.rcsb_pdbx_description,
            entity_poly_seq_length=r.entity_poly_seq_length,
        )
        for r in (ribo_struct.rnas or [])
    ]
    others = [
        PolymerEntry(
            auth_asym_id=o.auth_asym_id,
            assembly_id=o.assembly_id,
            nomenclature=_nomenclature_strings(o),
            rcsb_pdbx_description=o.rcsb_pdbx_description,
            entity_poly_seq_length=o.entity_poly_seq_length,
        )
        for o in (ribo_struct.other_polymers or [])
    ]

    assembly_map = None
    if ribo_struct.assembly_map:
        assembly_map = []
        for asm in ribo_struct.assembly_map:
            polys = [
                AssemblyPolymerInstance(
                    entity_id=inst.rcsb_polymer_entity_instance_container_identifiers.entity_id,
                    auth_asym_id=inst.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id,
                )
                for inst in (asm.polymer_entity_instances or [])
            ]
            nonpolys = None
            if asm.nonpolymer_entity_instances:
                nonpolys = [
                    AssemblyNonpolymerInstance(
                        entity_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.entity_id,
                        auth_asym_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.auth_asym_id,
                        auth_seq_id=inst.rcsb_nonpolymer_entity_instance_container_identifiers.auth_seq_id,
                    )
                    for inst in asm.nonpolymer_entity_instances
                ]
            assembly_map.append(AssemblyInfo(
                rcsb_id=asm.rcsb_id,
                polymer_entity_instances=polys,
                nonpolymer_entity_instances=nonpolys,
            ))

    return RibosomeProfile(
        rcsb_id=ribo_struct.rcsb_id,
        mitochondrial=bool(getattr(ribo_struct, "mitochondrial", False)),
        proteins=proteins,
        rnas=rnas,
        other_polymers=others,
        assembly_map=assembly_map,
    )


class RiboxyzStructureProvider:
    """Loads atoms + profile from local riboxyz assets (RibosomeOps)."""

    def fingerprint(self, rcsb_id: str) -> str:
        from ribctl.asset_manager.asset_types import AssetType
        p = AssetType.MMCIF.get_path(rcsb_id)
        return f"mmcif:{p}"

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        from ribctl.ribosome_ops import RibosomeOps
        from ribctl.asset_manager.asset_types import AssetType

        ro = RibosomeOps(rcsb_id)
        structure = ro.assets.biopython_structure()
        atoms = list(structure[0].get_atoms())
        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray([getattr(a, "element", "") or a.get_id()[0] for a in atoms])

        profile = _convert_profile(ro.profile)

        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(AssetType.MMCIF.get_path(rcsb_id)),
            "profile": profile,
            # Keep ro around for legacy stages that need biopython_structure etc.
            "ro": ro,
        }


class RiboxyzLandmarkProvider:
    def fingerprint(self, rcsb_id: str) -> str:
        return "ptc_via_trna+constriction_site:v1"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        from ribctl.lib.landmarks.ptc_via_trna import PTC_location
        from ribctl.lib.landmarks.constriction_site import get_constriction

        ptc = np.array(PTC_location(rcsb_id).location, dtype=np.float32)
        constr = np.array(get_constriction(rcsb_id), dtype=np.float32)
        return {"ptc_xyz": ptc, "constriction_xyz": constr}

```

ribctl/lib/npet2/adapters/standalone_providers.py
```py
# ribctl/lib/npet2/adapters/standalone_providers.py
"""
Standalone providers for running npet2 without the riboxyz repo.

Input sources:
  - mmCIF file on disk
  - Profile JSON (from riboxyz API or a local file)
  - Landmark coordinates (from riboxyz API or a local file)
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
from Bio.PDB.MMCIFParser import FastMMCIFParser

from ribctl.lib.npet2.core.ribosome_types import (
    ConstrictionInfo,
    PTCInfo,
    RibosomeProfile,
)
from ribctl.lib.npet2.core.settings import RIBOXYZ_API_BASE


class FileStructureProvider:
    """
    Load atoms from a local mmCIF file and profile from a JSON file or API.

    Usage:
        provider = FileStructureProvider(
            mmcif_path="/data/7K00.cif",
            profile_path="/data/7K00_profile.json",  # or None to fetch from API
        )
    """

    def __init__(
        self,
        mmcif_path: str | Path,
        profile_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.mmcif_path = Path(mmcif_path)
        self.profile_path = Path(profile_path) if profile_path else None
        self.api_base = api_base or RIBOXYZ_API_BASE

        if not self.mmcif_path.exists():
            raise FileNotFoundError(f"mmCIF file not found: {self.mmcif_path}")

    def fingerprint(self, rcsb_id: str) -> str:
        return f"mmcif:{self.mmcif_path}"

    def _load_profile(self, rcsb_id: str) -> RibosomeProfile:
        if self.profile_path and self.profile_path.exists():
            data = json.loads(self.profile_path.read_text())
            return RibosomeProfile.model_validate(data)

        # Fetch from API
        import requests

        url = f"{self.api_base}/structures/{rcsb_id.upper()}/profile"
        resp = requests.get(url, timeout=30)
        resp.raise_for_status()
        return RibosomeProfile.model_validate(resp.json())

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        parser = FastMMCIFParser(QUIET=True)
        structure = parser.get_structure(rcsb_id, str(self.mmcif_path))
        atoms = list(structure[0].get_atoms())

        if not atoms:
            raise ValueError(f"No atoms found in {self.mmcif_path}")

        xyz = np.asarray([a.get_coord() for a in atoms], dtype=np.float32)
        elem = np.asarray(
            [getattr(a, "element", "") or a.get_id()[0] for a in atoms]
        )

        profile = self._load_profile(rcsb_id)

        return {
            "atom_xyz": xyz,
            "atom_element": elem,
            "mmcif_path": str(self.mmcif_path),
            "profile": profile,
        }


class FileLandmarkProvider:
    """
    Load PTC/constriction from a local JSON file or the riboxyz API.

    Local file format:
        {
            "ptc": {"location": [x, y, z]},
            "constriction": {"location": [x, y, z]}
        }
    """

    def __init__(
        self,
        landmarks_path: Optional[str | Path] = None,
        api_base: Optional[str] = None,
    ):
        self.landmarks_path = Path(landmarks_path) if landmarks_path else None
        self.api_base = api_base or RIBOXYZ_API_BASE

    def fingerprint(self, rcsb_id: str) -> str:
        if self.landmarks_path:
            return f"landmarks_file:{self.landmarks_path}"
        return f"landmarks_api:{self.api_base}"

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        if self.landmarks_path and self.landmarks_path.exists():
            data = json.loads(self.landmarks_path.read_text())
            ptc_info = PTCInfo.model_validate(data["ptc"])
            constr_info = ConstrictionInfo.model_validate(data["constriction"])
            return {
                "ptc_xyz": np.array(ptc_info.location, dtype=np.float32),
                "constriction_xyz": np.array(constr_info.location, dtype=np.float32),
            }

        # Fetch from API
        import requests

        rcsb_id = rcsb_id.upper()

        ptc_resp = requests.get(
            f"{self.api_base}/loci/ptc", params={"rcsb_id": rcsb_id}, timeout=30
        )
        ptc_resp.raise_for_status()
        ptc_info = PTCInfo.model_validate(ptc_resp.json())

        constr_resp = requests.get(
            f"{self.api_base}/loci/constriction_site",
            params={"rcsb_id": rcsb_id},
            timeout=30,
        )
        constr_resp.raise_for_status()
        constr_info = ConstrictionInfo.model_validate(constr_resp.json())

        return {
            "ptc_xyz": np.array(ptc_info.location, dtype=np.float32),
            "constriction_xyz": np.array(constr_info.location, dtype=np.float32),
        }
```

ribctl/lib/npet2/backends/__init__.py
```py

```

ribctl/lib/npet2/backends/clustering_io.py
```py
from __future__ import annotations
from pathlib import Path
import json
import numpy as np

def clusters_from_labels(points: np.ndarray, labels: np.ndarray) -> dict[int, np.ndarray]:
    clusters: dict[int, list[int]] = {}
    for i, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(i)
    out = {}
    for lab, idxs in clusters.items():
        out[lab] = points[np.asarray(idxs, dtype=np.int32)]
    return out

def write_dbscan_pass(out_dir: Path, *, prefix: str, points: np.ndarray, labels: np.ndarray) -> dict:
    out_dir = out_dir / prefix
    out_dir.mkdir(parents=True, exist_ok=True)

    np.save(out_dir / "points.npy", points.astype(np.float32))
    np.save(out_dir / "labels.npy", labels.astype(np.int32))

    clusters = clusters_from_labels(points, labels)
    index = {"prefix": prefix, "n_points": int(points.shape[0]), "clusters": []}

    for cid, pts in sorted(clusters.items(), key=lambda kv: kv[0]):
        p = out_dir / f"cluster_id{cid}.npy"
        np.save(p, pts.astype(np.float32))
        index["clusters"].append({"id": int(cid), "n": int(pts.shape[0]), "path": p.name})

    (out_dir / "index.json").write_text(json.dumps(index, indent=2))
    return index

```

ribctl/lib/npet2/backends/grid_occupancy.py
```py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple

import numpy as np
from scipy import ndimage


@dataclass(frozen=True)
class GridSpec:
    # origin and voxel size define world coords: x = origin + i*voxel
    origin: np.ndarray          # (3,)
    voxel_size: float
    shape: Tuple[int, int, int] # (nx, ny, nz)


# ribctl/lib/npet2/backends/grid_occupancy.py

def make_cylinder_grid(radius_A: float, height_A: float, voxel_A: float,
                       z_min: float = 0.0) -> GridSpec:
    """
    Canonical cylinder in C0:
      x in [-R, R], y in [-R, R], z in [z_min, z_min + H]
    """
    nx = int(np.floor((2 * radius_A) / voxel_A)) + 1
    ny = int(np.floor((2 * radius_A) / voxel_A)) + 1
    nz = int(np.floor(height_A / voxel_A)) + 1
    origin = np.array([-radius_A, -radius_A, z_min], dtype=np.float32)
    return GridSpec(origin=origin, voxel_size=float(voxel_A), shape=(nx, ny, nz))


def grid_world_coords(grid: GridSpec) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Returns axis coordinate arrays (x, y, z) for voxel centers along each axis.
    """
    ox, oy, oz = grid.origin
    nx, ny, nz = grid.shape
    v = grid.voxel_size
    x = ox + np.arange(nx, dtype=np.float32) * v
    y = oy + np.arange(ny, dtype=np.float32) * v
    z = oz + np.arange(nz, dtype=np.float32) * v
    return x, y, z


def cylinder_mask(grid: GridSpec, radius_A: float) -> np.ndarray:
    """
    Boolean mask of voxels inside cylinder radius (in C0).
    """
    x, y, _ = grid_world_coords(grid)
    X, Y = np.meshgrid(x, y, indexing="ij")
    inside = (X * X + Y * Y) <= (radius_A * radius_A)
    # broadcast across z
    return inside[:, :, None]


def points_to_occupied_seeds(points_c0: np.ndarray, grid: GridSpec) -> np.ndarray:
    """
    Convert points (C0 coords) to a sparse occupied grid of seeds at nearest voxels.
    """
    pts = np.asarray(points_c0, dtype=np.float32)
    v = grid.voxel_size
    origin = grid.origin

    ijk = np.floor((pts - origin[None, :]) / v + 0.5).astype(np.int32)
    nx, ny, nz = grid.shape

    valid = (
        (ijk[:, 0] >= 0) & (ijk[:, 0] < nx) &
        (ijk[:, 1] >= 0) & (ijk[:, 1] < ny) &
        (ijk[:, 2] >= 0) & (ijk[:, 2] < nz)
    )
    ijk = ijk[valid]
    occ = np.zeros(grid.shape, dtype=np.bool_)
    if ijk.shape[0] > 0:
        occ[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True
    return occ


def occupancy_via_edt(points_c0: np.ndarray, grid: GridSpec, atom_radius_A: float) -> np.ndarray:
    """
    Occupancy grid: voxel is occupied if within atom_radius_A of any atom center.

    Steps:
      - seed occupied at nearest voxels
      - edt on ~occupied gives distance (in voxels) to nearest seed
      - threshold <= r_vox
    """
    seeds = points_to_occupied_seeds(points_c0, grid)

    # Distance (in voxels) from each voxel to nearest True in 'seeds'
    # distance_transform_edt computes distance to nearest zero;
    # so we compute on ~seeds, where zeros correspond to seeds.
    dist_vox = ndimage.distance_transform_edt(~seeds)

    r_vox = float(atom_radius_A) / float(grid.voxel_size)
    occupied = dist_vox <= r_vox
    return occupied


def empty_points_from_mask(grid: GridSpec, empty_mask: np.ndarray) -> np.ndarray:
    """
    Return coordinates (C0) of voxel centers where empty_mask is True.
    """
    empty_idx = np.where(empty_mask)
    x, y, z = grid_world_coords(grid)
    pts = np.column_stack((x[empty_idx[0]], y[empty_idx[1]], z[empty_idx[2]])).astype(np.float32)
    return pts


def save_grid_npy(grid: GridSpec, data: np.ndarray, path: Path, *, compress: bool = False) -> None:
    """
    Save a 3D grid along with its GridSpec metadata.
    File format: {path}_data.npy + {path}_spec.json
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if compress:
        np.savez_compressed(data_path.with_suffix('.npz'), data=data)
    else:
        np.save(data_path, data)
    
    spec_dict = {
        "origin": grid.origin.tolist(),
        "voxel_size": float(grid.voxel_size),
        "shape": list(grid.shape),
    }
    spec_path.write_text(__import__('json').dumps(spec_dict, indent=2))


def load_grid_npy(path: Path) -> tuple[GridSpec, np.ndarray]:
    """
    Load a grid saved by save_grid_npy.
    """
    data_path = path.parent / f"{path.stem}_data.npy"
    spec_path = path.parent / f"{path.stem}_spec.json"
    
    if data_path.with_suffix('.npz').exists():
        data = np.load(data_path.with_suffix('.npz'))['data']
    else:
        data = np.load(data_path)
    
    spec_dict = __import__('json').loads(spec_path.read_text())
    grid = GridSpec(
        origin=np.array(spec_dict["origin"], dtype=np.float32),
        voxel_size=float(spec_dict["voxel_size"]),
        shape=tuple(spec_dict["shape"]),
    )
    return grid, data


def voxel_to_world(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    """
    Convert voxel indices (i,j,k) to world coordinates.
    ijk: (N, 3) or (3,) array of voxel indices
    Returns: (N, 3) or (3,) world coordinates
    """
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin + ijk * grid.voxel_size


def world_to_voxel(grid: GridSpec, xyz: np.ndarray) -> np.ndarray:
    """
    Convert world coordinates to voxel indices.
    xyz: (N, 3) or (3,) world coordinates
    Returns: (N, 3) or (3,) voxel indices (floats; use floor/round as needed)
    """
    xyz = np.asarray(xyz, dtype=np.float32)
    return (xyz - grid.origin) / grid.voxel_size


def get_occupied_voxel_centers(grid: GridSpec, occupancy: np.ndarray) -> np.ndarray:
    """
    Get world coordinates of occupied voxel centers.
    """
    occupied_idx = np.argwhere(occupancy)
    return voxel_to_world(grid, occupied_idx)
# ribctl/lib/npet2/backends/grid_occupancy.py
# ADD these functions at the end:

from scipy import ndimage

def connected_components_3d(
    binary_mask: np.ndarray, 
    connectivity: int = 26
) -> tuple[np.ndarray, int]:
    """
    Find connected components in a 3D binary mask.
    
    Args:
        binary_mask: 3D boolean array
        connectivity: 6 (face), 18 (face+edge), or 26 (face+edge+corner)
    
    Returns:
        labeled: Array same shape as input with component labels (0=background)
        n_components: Number of components found
    """
    if connectivity == 6:
        structure = ndimage.generate_binary_structure(3, 1)
    elif connectivity == 18:
        structure = ndimage.generate_binary_structure(3, 2)
    elif connectivity == 26:
        structure = ndimage.generate_binary_structure(3, 3)
    else:
        raise ValueError(f"connectivity must be 6, 18, or 26, got {connectivity}")
    
    labeled, n_components = ndimage.label(binary_mask, structure=structure)
    return labeled, n_components


def get_largest_component(labeled: np.ndarray, n_components: int) -> np.ndarray:
    """
    Extract mask of the largest connected component.
    
    Args:
        labeled: Output from connected_components_3d
        n_components: Number of components
    
    Returns:
        Binary mask of largest component only
    """
    if n_components == 0:
        return np.zeros_like(labeled, dtype=bool)
    
    # Count voxels in each component (excluding background=0)
    component_sizes = np.bincount(labeled.ravel())
    component_sizes[0] = 0  # Ignore background
    
    largest_label = np.argmax(component_sizes)
    return labeled == largest_label


def get_component_stats(labeled: np.ndarray, n_components: int) -> list[dict]:
    """
    Get statistics for all connected components.
    
    Returns:
        List of dicts with {label, size, bbox_min, bbox_max}
    """
    stats = []
    
    for label in range(1, n_components + 1):
        mask = labeled == label
        size = int(mask.sum())
        
        if size == 0:
            continue
        
        indices = np.argwhere(mask)
        bbox_min = indices.min(axis=0)
        bbox_max = indices.max(axis=0)
        
        stats.append({
            "label": int(label),
            "size": size,
            "bbox_min": bbox_min.tolist(),
            "bbox_max": bbox_max.tolist(),
        })
    
    # Sort by size descending
    stats.sort(key=lambda x: x["size"], reverse=True)
    return stats


def morphological_clean(
    binary_mask: np.ndarray, 
    operation: str = "opening",
    iterations: int = 1
) -> np.ndarray:
    """
    Apply morphological operations to clean up a binary mask.
    
    Args:
        binary_mask: 3D boolean array
        operation: "opening" (remove small bits), "closing" (fill small holes), 
                  "erosion", "dilation"
        iterations: Number of times to apply operation
    
    Returns:
        Cleaned binary mask
    """
    if operation == "opening":
        return ndimage.binary_opening(binary_mask, iterations=iterations)
    elif operation == "closing":
        return ndimage.binary_closing(binary_mask, iterations=iterations)
    elif operation == "erosion":
        return ndimage.binary_erosion(binary_mask, iterations=iterations)
    elif operation == "dilation":
        return ndimage.binary_dilation(binary_mask, iterations=iterations)
    else:
        raise ValueError(f"Unknown operation: {operation}")
```

ribctl/lib/npet2/backends/meshing.py
```py
# ribctl/lib/npet2/backends/meshing.py
from __future__ import annotations
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.ndimage import binary_fill_holes, gaussian_filter

def save_mesh_with_ascii(mesh: pv.PolyData, path: Path, tag: str = "") -> None:
    """Save mesh as binary PLY + ASCII PLY side by side."""
    mesh.save(str(path))
    ascii_path = path.parent / f"{path.stem}_ascii.ply"
    try:
        mesh.save(str(ascii_path), binary=False)
    except Exception:
        try:
            import plyfile
            data = plyfile.PlyData.read(str(path))
            data.text = True
            data.write(str(ascii_path))
        except Exception as e:
            print(f"[meshing] ASCII PLY write failed{' (' + tag + ')' if tag else ''}: {e}")

def mesh_from_binary_volume(
    mask: np.ndarray,
    origin: np.ndarray,
    voxel_size: float,
    *,
    gaussian_sigma_voxels: float = 1.5,
    smooth_method: str = "taubin",
    smooth_iters: int = 20,
    taubin_pass_band: float = 0.1,
    fill_holes_size: float = 100.0,
) -> tuple[pv.PolyData, pv.PolyData]:
    """
    Marching cubes on a Gaussian-blurred binary volume, followed by mesh smoothing.

    Returns (smoothed_mesh, pre_smooth_mesh).
    Both are in the same coordinate frame as the input volume.
    Caller is responsible for coordinate transforms and saving.
    """
    vol = np.pad(mask.astype(np.float32), 2, constant_values=0.0)
    origin_pad = np.asarray(origin, dtype=np.float32) - 2 * voxel_size

    if gaussian_sigma_voxels > 0:
        vol = gaussian_filter(vol, sigma=gaussian_sigma_voxels)

    img = pv.ImageData(
        dimensions=vol.shape,
        spacing=(voxel_size, voxel_size, voxel_size),
        origin=(float(origin_pad[0]), float(origin_pad[1]), float(origin_pad[2])),
    )
    img.point_data["values"] = vol.ravel(order="F")

    surf = img.contour(isosurfaces=[0.5], scalars="values").triangulate()
    if surf.n_points == 0:
        raise ValueError("Marching cubes produced empty surface")

    surf = surf.clean(tolerance=0.0)

    if fill_holes_size > 0:
        surf = surf.fill_holes(fill_holes_size)

    surf = surf.connectivity(largest=True)

    pre_smooth = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)

    if smooth_iters > 0:
        if smooth_method == "taubin":
            surf = surf.smooth_taubin(n_iter=smooth_iters, pass_band=taubin_pass_band)
        else:
            surf = surf.smooth(n_iter=smooth_iters)

    surf = surf.compute_normals(auto_orient_normals=True, consistent_normals=True)
    return surf, pre_smooth


def voxelize_points(
    points: np.ndarray,
    voxel_size: float,
    pad_voxels: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    lo = points.min(axis=0) - pad_voxels * voxel_size
    hi = points.max(axis=0) + pad_voxels * voxel_size

    shape = tuple(np.ceil((hi - lo) / voxel_size).astype(int) + 1)

    ijk = np.floor((points - lo) / voxel_size + 0.5).astype(np.int32)
    for d in range(3):
        ijk[:, d] = np.clip(ijk[:, d], 0, shape[d] - 1)

    mask = np.zeros(shape, dtype=bool)
    mask[ijk[:, 0], ijk[:, 1], ijk[:, 2]] = True

    mask = binary_fill_holes(mask)
    return mask, lo.astype(np.float32)


def clip_mesh_to_atom_clearance(
    mesh: pv.PolyData,
    atom_xyz: np.ndarray,
    min_clearance_A: float = 1.5,
) -> pv.PolyData:
    from scipy.spatial import cKDTree

    tree = cKDTree(atom_xyz)
    pts = np.asarray(mesh.points, dtype=np.float64)

    dist, idx = tree.query(pts, k=1)
    violating = dist < min_clearance_A

    if violating.sum() == 0:
        return mesh

    nearest = atom_xyz[idx[violating]]
    direction = pts[violating] - nearest
    norms = np.maximum(np.linalg.norm(direction, axis=1, keepdims=True), 1e-8)
    pts[violating] = nearest + (direction / norms) * min_clearance_A

    result = mesh.copy()
    result.points = pts.astype(np.float32)
    return result
```

ribctl/lib/npet2/core/cache.py
```py
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import json
import shutil

from .run_id import stable_hash_dict

@dataclass(frozen=True)
class StageCacheKey:
    stage: str
    inputs_fp: dict
    params: dict
    impl_version: str = "v1.01"  # bump when you change semantics

    def digest(self) -> str:
        return stable_hash_dict({
            "stage": self.stage,
            "inputs_fp": self.inputs_fp,
            "params": self.params,
            "impl_version": self.impl_version,
        })[:20]

class LocalStageCache:
    def __init__(self, root: Path):
        self.root = root
        self.root.mkdir(parents=True, exist_ok=True)

    def entry_dir(self, key: StageCacheKey) -> Path:
        d = self.root / key.stage / key.digest()
        d.mkdir(parents=True, exist_ok=True)
        return d

    def has(self, key: StageCacheKey, required: list[str]) -> bool:
        d = self.root / key.stage / key.digest()
        if not d.exists():
            return False
        return all((d / r).exists() for r in required)

    def copy_into(self, key: StageCacheKey, dest: Path, files: list[str]) -> None:
        src = self.root / key.stage / key.digest()
        dest.mkdir(parents=True, exist_ok=True)
        for f in files:
            shutil.copy2(src / f, dest / f)

    def put_from(self, key: StageCacheKey, src_dir: Path, files: list[str]) -> None:
        dst = self.entry_dir(key)
        for f in files:
            shutil.copy2(src_dir / f, dst / f)

```

ribctl/lib/npet2/core/config.py
```py
# ribctl/lib/npet2/core/config.py
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Literal, Tuple


@dataclass(frozen=True)
class GridLevelConfig:

    name                 : str
    voxel_size_A         : float
    atom_radius_mode     : Literal["uniform", "vdw_bucket"] = "uniform"
    uniform_atom_radius_A: float = 2.0
    occupancy_backend    : Literal["legacy_kdtree", "edt"] = "legacy_kdtree"


@dataclass(frozen=True)
class RunConfig:
    # === Chain selection ===
    occupancy_chain_mode: Literal["walls_only", "assembly_all"] = "walls_only"
    occupancy_exclude_trna: bool = True
    occupancy_exclude_auth_asym_ids: Tuple[str, ...] = ()

    # === Region definition ===
    cylinder_radius_A: float = 35
    cylinder_height_A: float = 120
    cylinder_ptc_extension_A: float = 20

    # === Stage20: Exterior shell (whole ribosome surface) ===
    alpha_d3d_alpha       : float = 200
    alpha_d3d_tol         : float = 10
    alpha_d3d_offset      : float = 3
    alpha_kdtree_radius   : float = 40
    alpha_max_nn          : int   = 60
    alpha_tangent_planes_k: int   = 20
    alpha_poisson_depth   : int   = 6
    alpha_poisson_ptweight: int   = 4
    alpha_fill_holes      : float = 2000

    # === Stage40: Grid levels ===
    grid_levels: List[GridLevelConfig] = field(
        default_factory=lambda: [
            GridLevelConfig(
                name="level_0", voxel_size_A=1.0, occupancy_backend="legacy_kdtree"
            ),
        ]
    )

    # === Stage50: DBSCAN clustering on level_0 (1.0A grid) ===
    dbscan_level0_coarse_eps_A      : float = 5.5
    dbscan_level0_coarse_min_samples: int   = 600
    dbscan_level0_refine_eps_A      : float = 3.5
    dbscan_level0_refine_min_samples: int   = 175
    mesh_level0_enable              : bool  = True

    # === Stage55: Grid refinement (0.5A ROI pass) ===
    refine_voxel_size_A       : float = 0.5
    refine_roi_pad_A          : float = 10.0
    refine_atom_radius_A      : float = 2.0
    refine_keep_within_A      : float = 6.0
    refine_occ_close_iters    : int   = 0
    refine_void_open_iters    : int   = 1
    refine_forbid_roi_boundary: bool  = True

    # DBSCAN on refined grid
    dbscan_level1_coarse_eps_A      : float = 3.0
    dbscan_level1_coarse_min_samples: int   = 30
    dbscan_level1_refine_eps_A      : float = 3.0
    dbscan_level1_refine_min_samples: int   = 20

    refine_dbscan_max_points        : int   = 0
    refine_dbscan_seed              : int   = 0

    mesh_level1_enable: bool = True

    # === Meshing (MC + smoothing, shared by Stage50/55/70) ===
    mesh_smooth_method        : str   = "taubin"
    mesh_level0_gaussian_sigma: float = 1.0
    mesh_level1_gaussian_sigma: float = 1.5
    mesh_taubin_pass_band     : float = 0.1

    mesh_level0_smooth_iters  : int   = 40
    mesh_level1_smooth_iters  : int   = 60

    mesh_fill_holes_A         : float = 100.0
    mesh_atom_clearance_A     : float = 1.5
```

ribctl/lib/npet2/core/interfaces.py
```py
# ribctl/lib/npet2/core/interfaces.py
"""
Provider protocols for npet2.

These define the boundary between the pipeline and any data source.
Implement these to plug in riboxyz, a local file system, or an API.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Optional, Protocol

import numpy as np

from .ribosome_types import RibosomeProfile, PTCInfo, ConstrictionInfo
from .types import ArtifactRef, ArtifactType


class StructureProvider(Protocol):
    """Provides atom coordinates and the ribosome profile for a structure."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def load_atoms(self, rcsb_id: str) -> Dict[str, Any]:
        """
        Must return:
          - atom_xyz: (N, 3) float32
          - atom_element: (N,) str array (optional)
          - mmcif_path: str
          - profile: RibosomeProfile
        """
        ...


class LandmarkProvider(Protocol):
    """Provides PTC and constriction site coordinates."""

    def fingerprint(self, rcsb_id: str) -> str: ...

    def get_landmarks(self, rcsb_id: str) -> Dict[str, np.ndarray]:
        """
        Must return:
          - ptc_xyz: (3,) float32
          - constriction_xyz: (3,) float32
        """
        ...


class ArtifactStore(Protocol):
    @property
    def run_dir(self) -> Path: ...

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType,
                  data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_json(self, *, name: str, stage: str, obj: Any,
                 meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray,
                  meta: Optional[Dict[str, Any]] = None) -> ArtifactRef: ...

    def add_ref(self, ref: ArtifactRef) -> None: ...

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None: ...

```

ribctl/lib/npet2/core/manifest.py
```py
# ribctl/lib/npet2/core/manifest.py
from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import time

from .types import ArtifactRef, ArtifactType


@dataclass
class StageRecord:
    name: str
    status: str = "pending"  # pending|running|success|failure|skipped
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    params: Dict[str, Any] = field(default_factory=dict)
    note: Optional[str] = None


@dataclass
class RunManifest:
    rcsb_id: str
    run_id: str
    pipeline_version: str
    created_at: float = field(default_factory=lambda: time.time())

    inputs: Dict[str, Any] = field(default_factory=dict)
    config_resolved: Dict[str, Any] = field(default_factory=dict)

    stages: Dict[str, StageRecord] = field(default_factory=dict)
    artifacts: List[Dict[str, Any]] = field(default_factory=list)

    success: Optional[bool] = None
    error: Optional[str] = None

    def add_artifact(self, ref: ArtifactRef) -> None:
        self.artifacts.append({
            "name": ref.name,
            "type": ref.type.value,
            "path": str(ref.path),
            "stage": ref.stage,
            "meta": ref.meta,
            "depends_on": list(ref.depends_on),
        })

    def to_json(self) -> str:
        # dataclasses → dict
        d = asdict(self)
        # StageRecord needs manual flatten
        d["stages"] = {k: asdict(v) for k, v in self.stages.items()}
        return json.dumps(d, indent=2)

    @staticmethod
    def from_path(path: Path) -> "RunManifest":
        data = json.loads(path.read_text())
        m = RunManifest(
            rcsb_id=data["rcsb_id"],
            run_id=data["run_id"],
            pipeline_version=data["pipeline_version"],
            created_at=data.get("created_at", time.time()),
            inputs=data.get("inputs", {}),
            config_resolved=data.get("config_resolved", {}),
        )
        m.success = data.get("success")
        m.error = data.get("error")
        # stages
        for k, v in data.get("stages", {}).items():
            m.stages[k] = StageRecord(**v)
        m.artifacts = data.get("artifacts", [])
        return m

```

ribctl/lib/npet2/core/pipeline.py
```py
# ribctl/lib/npet2/core/pipeline.py
from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Dict, List
import time

from .types import StageContext


class Stage(ABC):
    key: str

    @abstractmethod
    def params(self, ctx: StageContext) -> Dict[str, Any]: ...

    @abstractmethod
    def run(self, ctx: StageContext) -> None: ...


def _fmt_value(v: Any) -> str:
    """Compact display for a parameter value."""
    if isinstance(v, float):
        # drop trailing zeros but keep one decimal
        return f"{v:g}"
    if isinstance(v, list) and len(v) > 3:
        return f"[{len(v)} items]"
    return str(v)


def _log_params(params: Dict[str, Any], prefix: str) -> None:
    if not params:
        return
    items = [f"{k}={_fmt_value(v)}" for k, v in params.items()
             if not isinstance(v, (dict, list))]
    # nested dicts/lists get their own lines
    nested = {k: v for k, v in params.items() if isinstance(v, (dict, list))}

    if items:
        line = ", ".join(items)
        # wrap at ~100 chars
        if len(line) > 100:
            mid = len(items) // 2
            print(f"  [{prefix}] {', '.join(items[:mid])}")
            print(f"  [{prefix}] {', '.join(items[mid:])}")
        else:
            print(f"  [{prefix}] {line}")
    for k, v in nested.items():
        if isinstance(v, list) and all(isinstance(x, dict) for x in v):
            for i, entry in enumerate(v):
                sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in entry.items())
                print(f"  [{prefix}] {k}[{i}]: {sub}")
        elif isinstance(v, dict):
            sub = ", ".join(f"{sk}={_fmt_value(sv)}" for sk, sv in v.items())
            print(f"  [{prefix}] {k}: {sub}")


class Pipeline:
    def __init__(self, stages: List[Stage]):
        self.stages = stages

    def run(self, ctx: StageContext) -> StageContext:
        n = len(self.stages)
        wall = 60

        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | run {ctx.run_id}")
        print(f"  stages: {n} | config: {type(ctx.config).__name__}")
        print("=" * wall)

        t_total = time.perf_counter()

        for i, stage in enumerate(self.stages, 1):
            params = stage.params(ctx)
            ctx.store.begin_stage(stage.key, params=params)

            print()
            print(f"--- [{i}/{n}] {stage.key} " + "-" * max(0, wall - len(stage.key) - 12))
            _log_params(params, stage.key)

            t0 = time.perf_counter()
            try:
                stage.run(ctx)
                dt = time.perf_counter() - t0
                print(f"  [{stage.key}] done in {dt:,.2f}s")
                ctx.store.end_stage(stage.key, success=True, note=f"elapsed_s={dt:.3f}")
            except Exception as e:
                dt = time.perf_counter() - t0
                print(f"  [{stage.key}] FAILED after {dt:,.2f}s: {e}")
                ctx.store.end_stage(
                    stage.key, success=False, note=f"elapsed_s={dt:.3f} err={e}"
                )
                ctx.store.finalize(success=False, error=str(e))
                raise

        dt_total = time.perf_counter() - t_total
        print()
        print("=" * wall)
        print(f"  npet2 | {ctx.rcsb_id} | completed in {dt_total:,.2f}s")
        print(f"  run_dir: {ctx.store.run_dir}")
        print("=" * wall)
        print()

        ctx.store.finalize(success=True)
        return ctx
```

ribctl/lib/npet2/core/polymer_enum.py
```py
# ribctl/lib/npet2/core/polymer_enum.py
"""
Ribosomal polymer nomenclature classes.

These enums encode the standard ribosomal nomenclature (Ban et al. 2014)
and are used by the pipeline for:
  - tRNA detection (to exclude from tunnel walls)
  - mL45 detection (mitochondrial tunnel debris)

They are intentionally duplicated from the outer riboxyz codebase so that
npet2 can run standalone. If you're integrating with riboxyz, the adapter
layer handles conversion.
"""
from __future__ import annotations

from enum import Enum
from typing import Union


class _PolymerEnumBase(str, Enum):
    def __repr__(self):
        return self.value


class tRNA(_PolymerEnumBase):
    tRNA = "tRNA"


class MitochondrialProteinClass(_PolymerEnumBase):
    # mSSU
    bS1m = "bS1m"; uS2m = "uS2m"; uS3m = "uS3m"; uS4m = "uS4m"
    uS5m = "uS5m"; bS6m = "bS6m"; uS7m = "uS7m"; uS8m = "uS8m"
    uS9m = "uS9m"; uS10m = "uS10m"; uS11m = "uS11m"; uS12m = "uS12m"
    uS13m = "uS13m"; uS14m = "uS14m"; uS15m = "uS15m"; bS16m = "bS16m"
    uS17m = "uS17m"; bS18m = "bS18m"; uS19m = "uS19m"; bS21m = "bS21m"
    mS22 = "mS22"; mS23 = "mS23"; mS25 = "mS25"; mS26 = "mS26"
    mS27 = "mS27"; mS29 = "mS29"; mS31 = "mS31"; mS33 = "mS33"
    mS34 = "mS34"; mS35 = "mS35"; mS37 = "mS37"; mS38 = "mS38"
    mS39 = "mS39"; mS40 = "mS40"; mS41 = "mS41"; mS42 = "mS42"
    mS43 = "mS43"; mS44 = "mS44"; mS45 = "mS45"; mS46 = "mS46"
    mS47 = "mS47"
    # mLSU
    uL1m = "uL1m"; uL2m = "uL2m"; uL3m = "uL3m"; uL4m = "uL4m"
    uL5m = "uL5m"; uL6m = "uL6m"; bL9m = "bL9m"; uL10m = "uL10m"
    uL11m = "uL11m"; bL12m = "bL12m"; uL13m = "uL13m"; uL14m = "uL14m"
    uL15m = "uL15m"; uL16m = "uL16m"; bL17m = "bL17m"; uL18m = "uL18m"
    bL19m = "bL19m"; bL20m = "bL20m"; bL21m = "bL21m"; uL22m = "uL22m"
    uL23m = "uL23m"; uL24m = "uL24m"; bL27m = "bL27m"; bL28m = "bL28m"
    uL29m = "uL29m"; uL30m = "uL30m"; bL31m = "bL31m"; bL32m = "bL32m"
    bL33m = "bL33m"; bL34m = "bL34m"; bL35m = "bL35m"; bL36m = "bL36m"
    mL37 = "mL37"; mL38 = "mL38"; mL39 = "mL39"; mL40 = "mL40"
    mL41 = "mL41"; mL42 = "mL42"; mL43 = "mL43"; mL44 = "mL44"
    mL45 = "mL45"; mL46 = "mL46"; mL48 = "mL48"; mL49 = "mL49"
    mL50 = "mL50"; mL51 = "mL51"; mL52 = "mL52"; mL53 = "mL53"
    mL54 = "mL54"; mL57 = "mL57"; mL58 = "mL58"; mL59 = "mL59"
    mL60 = "mL60"; mL61 = "mL61"; mL62 = "mL62"; mL63 = "mL63"
    mL64 = "mL64"; mL65 = "mL65"; mL66 = "mL66"; mL67 = "mL67"


class CytosolicProteinClass(_PolymerEnumBase):
    # SSU
    bS1 = "bS1"; eS1 = "eS1"; uS2 = "uS2"; uS3 = "uS3"
    uS4 = "uS4"; eS4 = "eS4"; uS5 = "uS5"; bS6 = "bS6"
    eS6 = "eS6"; uS7 = "uS7"; eS7 = "eS7"; uS8 = "uS8"
    eS8 = "eS8"; uS9 = "uS9"; uS10 = "uS10"; eS10 = "eS10"
    uS11 = "uS11"; uS12 = "uS12"; eS12 = "eS12"; uS13 = "uS13"
    uS14 = "uS14"; uS15 = "uS15"; bS16 = "bS16"; uS17 = "uS17"
    eS17 = "eS17"; bS18 = "bS18"; uS19 = "uS19"; eS19 = "eS19"
    bS20 = "bS20"; bS21 = "bS21"; bTHX = "bTHX"; eS21 = "eS21"
    eS24 = "eS24"; eS25 = "eS25"; eS26 = "eS26"; eS27 = "eS27"
    eS28 = "eS28"; eS30 = "eS30"; eS31 = "eS31"; RACK1 = "RACK1"
    # LSU
    uL1 = "uL1"; uL2 = "uL2"; uL3 = "uL3"; uL4 = "uL4"
    uL5 = "uL5"; uL6 = "uL6"; eL6 = "eL6"; eL8 = "eL8"
    bL9 = "bL9"; uL10 = "uL10"; uL11 = "uL11"; bL12 = "bL12"
    uL13 = "uL13"; eL13 = "eL13"; uL14 = "uL14"; eL14 = "eL14"
    uL15 = "uL15"; eL15 = "eL15"; uL16 = "uL16"; bL17 = "bL17"
    uL18 = "uL18"; eL18 = "eL18"; bL19 = "bL19"; eL19 = "eL19"
    bL20 = "bL20"; eL20 = "eL20"; bL21 = "bL21"; eL21 = "eL21"
    uL22 = "uL22"; eL22 = "eL22"; uL23 = "uL23"; uL24 = "uL24"
    eL24 = "eL24"; bL25 = "bL25"; bL27 = "bL27"; eL27 = "eL27"
    bL28 = "bL28"; eL28 = "eL28"; uL29 = "uL29"; eL29 = "eL29"
    uL30 = "uL30"; eL30 = "eL30"; bL31 = "bL31"; eL31 = "eL31"
    bL32 = "bL32"; eL32 = "eL32"; bL33 = "bL33"; eL33 = "eL33"
    bL34 = "bL34"; eL34 = "eL34"; bL35 = "bL35"; bL36 = "bL36"
    eL36 = "eL36"; eL37 = "eL37"; eL38 = "eL38"; eL39 = "eL39"
    eL40 = "eL40"; eL41 = "eL41"; eL42 = "eL42"; eL43 = "eL43"
    P1P2 = "P1P2"


class MitochondrialRNAClass(_PolymerEnumBase):
    mtrRNA12S = "mt12SrRNA"
    mtrRNA16S = "mt16SrRNA"


class CytosolicRNAClass(_PolymerEnumBase):
    rRNA_5S   = "5SrRNA"
    rRNA_16S  = "16SrRNA"
    rRNA_23S  = "23SrRNA"
    rRNA_25S  = "25SrRNA"
    rRNA_5_8S = "5.8SrRNA"
    rRNA_18S  = "18SrRNA"
    rRNA_28S  = "28SrRNA"


class ElongationFactorClass(_PolymerEnumBase):
    eEF1A = "eEF1A"; eEF1B = "eEF1B"; eFSec = "eFSec"; eEF2 = "eEF2"
    mtEF4 = "mtEF4"; eIF5A = "eIF5A"; eEF3 = "eEF3"
    EF_Tu = "EF-Tu"; EF_Ts = "EF-Ts"; SelB = "SelB"; EF_G = "EF-G"
    EF4 = "EF4"; EF_P = "EF-P"; Tet_O = "Tet_O"; Tet_M = "Tet_M"
    RelA = "RelA"; BipA = "BipA"
    aEF1A = "aEF1A"; aEF2 = "aEF2"


class InitiationFactorClass(_PolymerEnumBase):
    eIF1 = "eIF1"; eIF1A = "eIF1A"
    eIF2_alpha = "eIF2_alpha"; eIF2_beta = "eIF2_beta"; eIF2_gamma = "eIF2_gamma"
    eIF2B_alpha = "eIF2B_alpha"; eIF2B_beta = "eIF2B_beta"
    eIF2B_gamma = "eIF2B_gamma"; eIF2B_delta = "eIF2B_delta"; eIF2B_epsilon = "eIF2B_epsilon"
    eIF3_subunitA = "eIF3_subunitA"; eIF3_subunitB = "eIF3_subunitB"
    eIF3_subunitC = "eIF3_subunitC"; eIF3_subunitD = "eIF3_subunitD"
    eIF3_subunitE = "eIF3_subunitE"; eIF3_subunitF = "eIF3_subunitF"
    eIF3_subunitG = "eIF3_subunitG"; eIF3_subunitH = "eIF3_subunitH"
    eIF3_subunitI = "eIF3_subunitI"; eIF3_subunitJ = "eIF3_subunitJ"
    eIF3_subunitK = "eIF3_subunitK"; eIF3_subunitL = "eIF3_subunitL"
    eIF3_subunitM = "eIF3_subunitM"
    eIF4F_4A = "eIF4F_4A"; eIF4F_4G = "eIF4F_4G"; eIF4F_4E = "eIF4F_4E"
    eIF4B = "eIF4B"; eIF5B = "eIF5B"; eIF5 = "eIF5"
    IF1 = "IF1"; IF2 = "IF2"; IF3 = "IF3"
    aIF_1A = "aIF1A"; aIF_2_alpha = "aIF2_alpha"; aIF_2_beta = "aIF2_beta"
    aIF_2_gamma = "aIF2_gamma"; aIF_2B_alpha = "aIF2B_alpha"
    aIF_2B_beta = "aIF2B_beta"; aIF_2B_delta = "aIF2B_delta"
    aIF5A = "aIF5A"; aIF5B = "aIF5B"


# Composite types
ProteinClass = Union[MitochondrialProteinClass, CytosolicProteinClass]
LifecycleFactorClass = Union[ElongationFactorClass, InitiationFactorClass]
PolypeptideClass = Union[LifecycleFactorClass, ProteinClass]
PolynucleotideClass = Union[CytosolicRNAClass, MitochondrialRNAClass, tRNA]
PolymerClass = Union[PolynucleotideClass, PolypeptideClass]


def parse_polymer_class(value: str) -> PolymerClass:
    """Parse a string into the appropriate PolymerClass enum member.

    Used when deserializing profile JSON where nomenclature entries are strings.
    Raises ValueError if the string doesn't match any known polymer class.
    """
    for enum_cls in (
        tRNA,
        CytosolicProteinClass, MitochondrialProteinClass,
        CytosolicRNAClass, MitochondrialRNAClass,
        ElongationFactorClass, InitiationFactorClass,
    ):
        try:
            return enum_cls(value)
        except ValueError:
            continue
    raise ValueError(f"Unknown polymer class: {value!r}")
```

ribctl/lib/npet2/core/ribosome_types.py
```py
# ribctl/lib/npet2/core/ribosome_types.py
"""
Standalone pydantic models for ribosome structure profiles as consumed by npet2.

These mirror the relevant subset of the riboxyz schema. When running inside the
riboxyz repo, the adapter converts full RibosomeStructure objects to these.
When running standalone, users provide profile JSON that validates against these
models directly.
"""
from __future__ import annotations

from typing import Optional

from pydantic import BaseModel, field_validator

from .polymer_enum import PolymerClass, parse_polymer_class


# ---------------------------------------------------------------------------
# Polymer entries (slim -- only what the pipeline touches)
# ---------------------------------------------------------------------------

class PolymerEntry(BaseModel):
    """A single polymer chain as seen by the npet2 pipeline."""
    auth_asym_id: str
    assembly_id: int = 0
    nomenclature: list[PolymerClass] = []
    rcsb_pdbx_description: Optional[str] = None
    entity_poly_seq_length: int = 0

    @field_validator("nomenclature", mode="before")
    @classmethod
    def _coerce_nomenclature(cls, v):
        """Accept raw strings (from JSON) and convert to enum members."""
        if not v:
            return []
        out = []
        for item in v:
            if isinstance(item, str):
                out.append(parse_polymer_class(item))
            else:
                out.append(item)
        return out


class ProteinEntry(PolymerEntry):
    """Protein chain. Inherits all pipeline-relevant fields from PolymerEntry."""
    pass


class RNAEntry(PolymerEntry):
    """RNA chain."""
    pass


# ---------------------------------------------------------------------------
# Assembly map (for first-assembly filtering)
# ---------------------------------------------------------------------------

class AssemblyPolymerInstance(BaseModel):
    entity_id: str
    auth_asym_id: str


class AssemblyNonpolymerInstance(BaseModel):
    entity_id: str
    auth_asym_id: str
    auth_seq_id: str = ""


class AssemblyInfo(BaseModel):
    rcsb_id: str
    polymer_entity_instances: list[AssemblyPolymerInstance] = []
    nonpolymer_entity_instances: Optional[list[AssemblyNonpolymerInstance]] = None


# ---------------------------------------------------------------------------
# Top-level profile
# ---------------------------------------------------------------------------

class RibosomeProfile(BaseModel):
    """
    Minimal ribosome structure profile for the npet2 pipeline.

    This is the only "ribosome knowledge" npet2 needs beyond the mmCIF itself
    and the landmark coordinates. It tells the pipeline which chains form the
    tunnel walls (proteins + rRNAs), which to exclude (tRNAs, factors), and
    whether the structure is mitochondrial.

    Can be loaded from:
      - riboxyz API: GET /structures/{rcsb_id}/profile
      - Local JSON file (--profile path/to/profile.json)
      - Constructed programmatically by the riboxyz adapter
    """
    rcsb_id: str
    mitochondrial: bool = False

    proteins: list[ProteinEntry] = []
    rnas: list[RNAEntry] = []
    other_polymers: list[PolymerEntry] = []

    assembly_map: Optional[list[AssemblyInfo]] = None

    def all_polymers(self) -> list[PolymerEntry]:
        return [*self.proteins, *self.rnas, *self.other_polymers]

    def first_assembly_auth_asym_ids(self) -> list[str]:
        """Return auth_asym_ids from the first assembly, or raise."""
        if not self.assembly_map:
            raise ValueError(
                f"No assembly_map in profile for {self.rcsb_id}. "
                "Provide one, or set occupancy_chain_mode='walls_only' to skip assembly filtering."
            )
        first = self.assembly_map[0]
        ids = [inst.auth_asym_id for inst in first.polymer_entity_instances]
        if first.nonpolymer_entity_instances:
            ids.extend(inst.auth_asym_id for inst in first.nonpolymer_entity_instances)
        return ids


# ---------------------------------------------------------------------------
# Landmark types
# ---------------------------------------------------------------------------

class PTCInfo(BaseModel):
    """Peptidyl transferase center location."""
    location: list[float]

    @field_validator("location")
    @classmethod
    def _validate_location(cls, v):
        if len(v) != 3:
            raise ValueError(f"PTC location must be [x, y, z], got {len(v)} values")
        return [float(x) for x in v]


class ConstrictionInfo(BaseModel):
    """Constriction site location."""
    location: list[float]

    @field_validator("location")
    @classmethod
    def _validate_location(cls, v):
        if len(v) != 3:
            raise ValueError(f"Constriction location must be [x, y, z], got {len(v)} values")
        return [float(x) for x in v]
```

ribctl/lib/npet2/core/run_id.py
```py
# ribctl/lib/npet2/core/run_id.py
from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime
from pathlib import Path
from typing import Any, Dict


def stable_hash_dict(d: Dict[str, Any]) -> str:
    payload = json.dumps(d, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _next_seq_index(runs_dir: Path) -> int:
    """
    Scan existing run directories under runs_dir for the pattern NNN_...
    and return max+1.  If none exist, returns 1.
    """
    max_idx = 0
    if runs_dir.exists():
        for d in runs_dir.iterdir():
            if d.is_dir():
                m = re.match(r"^(\d{3,})_", d.name)
                if m:
                    max_idx = max(max_idx, int(m.group(1)))
    return max_idx + 1


def compute_run_id(
    *,
    rcsb_id: str,
    pipeline_version: str,
    inputs_fp: Dict[str, str],
    config_resolved: Dict[str, Any],
    runs_dir: Path,
) -> str:
    """
    run_id = SEQ_TIMESTAMP_HASH

    Format: NNN_YYYYMMDD_HHMMSS_<hash16>
    Sequential index makes it trivial to find latest run.
    """
    blob = {
        "rcsb_id": rcsb_id.upper(),
        "pipeline_version": pipeline_version,
        "inputs": dict(sorted(inputs_fp.items())),
        "config": config_resolved,
    }
    hash_str = stable_hash_dict(blob)[:16]
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    seq = _next_seq_index(runs_dir)
    return f"{seq:03d}_{timestamp}_{hash_str}"

```

ribctl/lib/npet2/core/settings.py
```py
# ribctl/lib/npet2/core/settings.py
"""
Default paths and external tool locations for npet2.

All of these can be overridden via:
  - Environment variables (NPET2_ROOT, NPET2_POISSON_RECON_BIN, etc.)
  - CLI flags
  - RunConfig / programmatic construction
"""
from __future__ import annotations

import os
from pathlib import Path


def _env_path(var: str, default: str) -> Path:
    return Path(os.environ.get(var, default))


NPET2_ROOT        = _env_path("NPET2_ROOT", str(Path.home() / "npet2_data"))
NPET2_RUNS_ROOT   = _env_path("NPET2_RUNS_ROOT", str(NPET2_ROOT / "runs"))
NPET2_CACHE_ROOT  = _env_path("NPET2_CACHE_ROOT", str(NPET2_ROOT / "cache"))
POISSON_RECON_BIN = os.environ.get("NPET2_POISSON_RECON_BIN", "PoissonRecon")
RIBOXYZ_API_BASE  = os.environ.get("NPET2_RIBOXYZ_API_URL", "http://localhost:8000")

```

ribctl/lib/npet2/core/store.py
```py
# ribctl/lib/npet2/core/store.py
from __future__ import annotations

import json
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

from .interfaces import ArtifactStore
from .manifest import RunManifest, StageRecord
from .types import ArtifactRef, ArtifactType


class LocalRunStore(ArtifactStore):
    def __init__(self, run_dir: Path, manifest: RunManifest):
        self._run_dir       = run_dir
        self._manifest      = manifest
        self._manifest_path = run_dir / "manifest.json"
        self._run_dir.mkdir(parents=True, exist_ok=True)
        self._write_manifest()

    def _abs(self, p: Path) -> Path:
        return p if p.is_absolute() else (self.run_dir / p)

    def _rel(self, p: Path) -> str:
        p = self._abs(p)
        try:
            return str(p.relative_to(self.run_dir))
        except ValueError:
            # Not under run_dir; fall back to absolute string (still tracked)
            return str(p)


    def register_file(
        self,
        *,
        name: str,
        stage: str,
        type: ArtifactType,
        path: Path,
        meta: Optional[Dict[str, Any]] = None,
        depends_on: tuple[str, ...] = (),
    ) -> ArtifactRef:
        ap = self._abs(path)

        ref = ArtifactRef(
            name=name,
            type=type,
            path=self._rel(ap),
            stage=stage,
            meta=meta or {},
            depends_on=depends_on,
        )
        self.add_ref(ref)
        return ref


    @property
    def run_dir(self) -> Path:
        return self._run_dir

    @property
    def manifest(self) -> RunManifest:
        return self._manifest

    def stage_dir(self, stage: str) -> Path:
        d = self._run_dir / "stage" / stage
        d.mkdir(parents=True, exist_ok=True)
        return d

    def _write_manifest(self) -> None:
        self._manifest_path.write_text(self._manifest.to_json())

    def add_ref(self, ref: ArtifactRef) -> None:
        self._manifest.add_artifact(ref)
        self._write_manifest()

    def put_bytes(self, *, name: str, stage: str, type: ArtifactType, data: bytes, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        suffix = {
            ArtifactType.JSON: ".json",
            ArtifactType.NUMPY: ".npy",
            ArtifactType.PNG: ".png",
            ArtifactType.TXT: ".txt",
            ArtifactType.PLY_MESH: ".ply",
            ArtifactType.PLY_PCD: ".ply",
        }[type]
        out = self.stage_dir(stage) / f"{name}{suffix}"
        out.write_bytes(data)
        ref = ArtifactRef(name=name, type=type, path=out, stage=stage, meta=meta or {})
        self.add_ref(ref)
        return ref

    def put_json(self, *, name: str, stage: str, obj: Any, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        data = json.dumps(obj, indent=2).encode("utf-8")
        return self.put_bytes(name=name, stage=stage, type=ArtifactType.JSON, data=data, meta=meta)

    def put_numpy(self, *, name: str, stage: str, arr: np.ndarray, meta: Optional[Dict[str, Any]] = None) -> ArtifactRef:
        out = self.stage_dir(stage) / f"{name}.npy"
        np.save(out, arr)
        ref = ArtifactRef(name=name, type=ArtifactType.NUMPY, path=out, stage=stage, meta=meta or {})
        self.add_ref(ref)
        return ref

    # Stage status helpers (optional but useful)
    def begin_stage(self, stage: str, params: Optional[Dict[str, Any]] = None) -> None:
        rec = self._manifest.stages.get(stage) or StageRecord(name=stage)
        rec.status = "running"
        rec.started_at = time.time()
        rec.params = params or {}
        self._manifest.stages[stage] = rec
        self._write_manifest()

    def end_stage(self, stage: str, success: bool, note: Optional[str] = None) -> None:
        rec = self._manifest.stages.get(stage) or StageRecord(name=stage)
        rec.status = "success" if success else "failure"
        rec.ended_at = time.time()
        rec.note = note
        self._manifest.stages[stage] = rec
        self._write_manifest()

    def finalize(self, *, success: bool, error: Optional[str] = None) -> None:
        self._manifest.success = success
        self._manifest.error = error
        self._write_manifest()

```

ribctl/lib/npet2/core/structure_selection.py
```py
# ribctl/lib/npet2/core/structure_selection.py
"""
Chain selection policies for tunnel wall definition.

Works with npet2-internal RibosomeProfile types exclusively.
"""

from __future__ import annotations

from typing import Iterable, Set

from .ribosome_types import RibosomeProfile, PolymerEntry
from .polymer_enum import tRNA as tRNAClass

_TRNA_HINTS = ("trna", "transfer rna")


def looks_like_trna(poly: PolymerEntry) -> bool:
    """Heuristic tRNA detection from nomenclature and description."""
    for nom in poly.nomenclature:
        if isinstance(nom, tRNAClass):
            return True
        if "trna" in str(nom).lower():
            return True

    desc = (poly.rcsb_pdbx_description or "").lower()
    if any(h in desc for h in _TRNA_HINTS):
        return True

    return False


def ribosome_wall_auth_asym_ids(
    profile: RibosomeProfile,
    *,
    exclude_trna: bool = True,
    extra_exclude: Iterable[str] = (),
) -> Set[str]:
    """
    Tunnel wall = ribosomal proteins + rRNAs.

    This automatically excludes waters/ions/ligands/nonpolymers (they aren't
    in proteins or rnas). Modified residues within ribosomal polymers are
    included because they are covalently part of the wall.
    """
    wall = {p.auth_asym_id for p in profile.proteins}
    wall |= {r.auth_asym_id for r in profile.rnas}

    if exclude_trna:
        all_polys = profile.all_polymers()
        trna_ids = {p.auth_asym_id for p in all_polys if looks_like_trna(p)}
        wall -= trna_ids

    wall -= set(extra_exclude)
    return wall


def intersect_with_first_assembly(
    profile: RibosomeProfile, chain_ids: Set[str]
) -> Set[str]:
    """Only keep chains present in the first assembly (if assembly_map available)."""
    try:
        asm_ids = set(profile.first_assembly_auth_asym_ids())
        return chain_ids & asm_ids
    except (ValueError, IndexError):
        return chain_ids


def tunnel_debris_chains(rcsb_id: str, profile: RibosomeProfile) -> list[str]:
    """
    Hardcoded per-structure chain exclusions for known tunnel debris.

    Also handles mitochondrial mL45 (sits inside the exit tunnel).
    """
    from .polymer_enum import MitochondrialProteinClass

    DEBRIS_MAP = {
        "3J7Z": ["a", "7"],
        "5GAK": ["z"],
        "5NWY": ["s"],
        "7A5G": ["Y2"],
        "9F1D": ["BK"],
    }
    skip = list(DEBRIS_MAP.get(rcsb_id.upper(), []))

    if profile.mitochondrial:
        for poly in profile.all_polymers():
            if MitochondrialProteinClass.mL45 in poly.nomenclature:
                skip.append(poly.auth_asym_id)
                break

    return skip


def atom_inclusion_policy(
    profile: RibosomeProfile,
    config,
    rcsb_id: str,
) -> dict:
    """
    Central policy for which atoms go into occupancy calculations.

    Returns:
        wall_chain_ids: set of auth_asym_ids for occupancy
        excluded_chain_ids: set of auth_asym_ids removed
        reasons: dict mapping excluded chain_id -> reason string
    """
    debris = tunnel_debris_chains(rcsb_id, profile)
    manual_exclude = list(getattr(config, "occupancy_exclude_auth_asym_ids", ()))
    exclude_trna = bool(getattr(config, "occupancy_exclude_trna", True))

    all_exclude = list(dict.fromkeys(debris + manual_exclude))

    wall = ribosome_wall_auth_asym_ids(
        profile, exclude_trna=exclude_trna, extra_exclude=all_exclude
    )
    wall = intersect_with_first_assembly(profile, wall)

    reasons = {}
    for c in debris:
        reasons[c] = "tunnel_debris (hardcoded)"
    for c in manual_exclude:
        if c not in reasons:
            reasons[c] = "config exclude"
    if exclude_trna:
        for p in profile.all_polymers():
            if looks_like_trna(p) and p.auth_asym_id not in wall:
                reasons[p.auth_asym_id] = "tRNA (auto-detected)"

    return {
        "wall_chain_ids": wall,
        "excluded_chain_ids": set(reasons.keys()),
        "reasons": reasons,
    }

```

ribctl/lib/npet2/core/types.py
```py
# ribctl/lib/npet2/core/types.py
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Mapping, Optional


class ArtifactType(str, Enum):
    JSON = "json"
    NUMPY = "npy"
    PLY_MESH = "ply_mesh"
    PLY_PCD = "ply_pcd"
    PNG = "png"
    TXT = "txt"


@dataclass(frozen=True)
class ArtifactRef:
    """
    A stable handle to an on-disk artifact, referenced from the manifest.
    """
    name: str                 # semantic name: "ptc", "empty_points_level_0"
    type: ArtifactType
    path: Path                # absolute or run-relative; store decides
    stage: str                # stage key: "10_landmarks"
    meta: Dict[str, Any] = field(default_factory=dict)
    depends_on: tuple[str, ...] = ()  # artifact names (or ids later)


@dataclass
class StageContext:
    """
    Shared context passed through the pipeline.
    - inputs: raw objects needed by compute (arrays, coords, providers, etc.)
    - artifacts: ArtifactRef registry for cross-stage access
    - stats: cheap summaries to help decisions (counts/bounds, etc.)
    """
    run_id: str
    rcsb_id: str
    config: Any  # RunConfig (kept Any to avoid import cycles)
    store: Any   # ArtifactStore

    inputs: Dict[str, Any] = field(default_factory=dict)
    artifacts: Dict[str, ArtifactRef] = field(default_factory=dict)
    stats: Dict[str, Any] = field(default_factory=dict)

    def require(self, key: str) -> Any:
        if key not in self.inputs:
            raise KeyError(f"Missing required input: {key}")
        return self.inputs[key]

    def require_artifact(self, name: str) -> ArtifactRef:
        if name not in self.artifacts:
            raise KeyError(f"Missing required artifact: {name}")
        return self.artifacts[name]

```

ribctl/lib/npet2/stages/bootstrap.py
```py
# ribctl/lib/npet2/stages/bootstrap.py (updated)
from __future__ import annotations

from typing import Any, Dict

import numpy as np

from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.ribosome_types import RibosomeProfile
from ribctl.lib.npet2.core.types import StageContext


class Stage00Inputs(Stage):
    key = "00_inputs"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        structure_provider = ctx.require("structure_provider")
        data = structure_provider.load_atoms(ctx.rcsb_id)

        atom_xyz = np.asarray(data["atom_xyz"], dtype=np.float32)
        ctx.inputs["atom_xyz"] = atom_xyz
        ctx.inputs["atom_element"] = data.get("atom_element", None)
        ctx.inputs["mmcif_path"] = data["mmcif_path"]

        # Profile: validate it's the right type
        profile = data["profile"]
        if not isinstance(profile, RibosomeProfile):
            profile = RibosomeProfile.model_validate(profile)
        ctx.inputs["profile"] = profile

        # If the provider gave us a biopython structure or RibosomeOps, stash it.
        # Standalone providers won't -- Stage30 will parse mmcif on demand.
        if "ro" in data:
            ro = data["ro"]
            ctx.inputs["ro"] = ro
            ctx.inputs["biopython_structure"] = ro.assets.biopython_structure()
        elif "biopython_structure" in data:
            ctx.inputs["biopython_structure"] = data["biopython_structure"]

        ctx.artifacts["atom_xyz"] = ctx.store.put_numpy(
            name="atom_xyz",
            stage=self.key,
            arr=atom_xyz,
            meta={"shape": list(atom_xyz.shape), "dtype": str(atom_xyz.dtype)},
        )

        mins = atom_xyz.min(axis=0)
        maxs = atom_xyz.max(axis=0)
        ctx.stats["atom_bounds"] = {"min": mins.tolist(), "max": maxs.tolist()}
        ctx.stats["n_atoms"] = int(atom_xyz.shape[0])


class Stage10Landmarks(Stage):
    key = "10_landmarks"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        landmark_provider = ctx.require("landmark_provider")
        lm = landmark_provider.get_landmarks(ctx.rcsb_id)

        ptc = np.asarray(lm["ptc_xyz"], dtype=np.float32)
        constr = np.asarray(lm["constriction_xyz"], dtype=np.float32)

        if ptc.shape != (3,):
            raise ValueError(f"PTC must be shape (3,), got {ptc.shape}")
        if constr.shape != (3,):
            raise ValueError(f"Constriction must be shape (3,), got {constr.shape}")

        ctx.inputs["ptc_xyz"] = ptc
        ctx.inputs["constriction_xyz"] = constr

        D = float(np.linalg.norm(constr - ptc))
        if D < 1.0:
            raise ValueError(
                f"PTC and constriction are too close ({D:.1f}A) -- check coordinates"
            )

        z_min = -float(ctx.config.cylinder_ptc_extension_A)
        z_max = float(ctx.config.cylinder_height_A)

        ctx.inputs["cylinder_z_min"] = z_min
        ctx.inputs["cylinder_z_max"] = z_max
        ctx.inputs["landmark_distance"] = D

        print(
            f"  [10_landmarks] PTC-Constriction distance={D:.1f}A, "
            f"cylinder z=[{z_min:.1f}, {z_max:.1f}]A"
        )

        ctx.artifacts["ptc"] = ctx.store.put_json(
            name="ptc",
            stage=self.key,
            obj={"location": ptc.tolist()},
            meta={"units": "A"},
        )
        ctx.artifacts["constriction_site"] = ctx.store.put_json(
            name="constriction_site",
            stage=self.key,
            obj={"location": constr.tolist()},
            meta={
                "units": "A",
                "landmark_distance_A": D,
                "cylinder_z_min": z_min,
                "cylinder_z_max": z_max,
            },
        )

```

ribctl/lib/npet2/stages/grid_refine.py
```py
# ribctl/lib/npet2/stages/grid_refine.py

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from typing import Any, Dict, Tuple

import numpy as np
import pyvista as pv
from scipy import ndimage
import time
import open3d as o3d

from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from ribctl.lib.npet.kdtree_approach import (
    transform_points_to_C0,
    transform_points_from_C0,
    estimate_normals,
)

from ribctl.lib.npet2.backends.grid_occupancy import (
    GridSpec,
    occupancy_via_edt,
)


def _make_bbox_grid(lo: np.ndarray, hi: np.ndarray, voxel: float) -> GridSpec:
    lo = np.asarray(lo, dtype=np.float32)
    hi = np.asarray(hi, dtype=np.float32)
    voxel = float(voxel)

    span = hi - lo
    shape = tuple((np.ceil(span / voxel).astype(np.int32) + 1).tolist())
    return GridSpec(origin=lo, voxel_size=voxel, shape=shape)


def _voxel_centers_from_indices(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    ijk = np.asarray(ijk, dtype=np.float32)
    return grid.origin[None, :] + ijk * float(grid.voxel_size)


def _points_to_ijk(grid: GridSpec, pts_c0: np.ndarray) -> np.ndarray:
    """Nearest-voxel mapping for points in C0 -> ijk indices."""
    v = float(grid.voxel_size)
    ijk = np.floor((pts_c0 - grid.origin[None, :]) / v + 0.5).astype(np.int32)
    return ijk


def _valid_ijk(grid: GridSpec, ijk: np.ndarray) -> np.ndarray:
    nx, ny, nz = grid.shape
    m = (
        (ijk[:, 0] >= 0)
        & (ijk[:, 0] < nx)
        & (ijk[:, 1] >= 0)
        & (ijk[:, 1] < ny)
        & (ijk[:, 2] >= 0)
        & (ijk[:, 2] < nz)
    )
    return m




class Stage55GridRefine(Stage):
    key = "55_grid_refine"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "voxel_size_A"         : float(c.refine_voxel_size_A),
            "roi_pad_A"            : float(c.refine_roi_pad_A),
            "atom_radius_A"        : float(c.refine_atom_radius_A),
            "keep_within_A"        : float(c.refine_keep_within_A),
            "occ_close_iters"      : int(c.refine_occ_close_iters),
            "void_open_iters"      : int(c.refine_void_open_iters),
            "forbid_roi_boundary"  : bool(c.refine_forbid_roi_boundary),
            "coarse_eps_A"         : float(c.dbscan_level1_coarse_eps_A),
            "coarse_min_samples"   : int(c.dbscan_level1_coarse_min_samples),
            "refine_eps_A"         : float(c.dbscan_level1_refine_eps_A),
            "refine_min_samples"   : int(c.dbscan_level1_refine_min_samples),
            "dbscan_max_points"    : int(c.refine_dbscan_max_points),
            "dbscan_seed"          : int(c.refine_dbscan_seed),
            "mesh_enable"          : bool(getattr(c, "mesh_level1_enable", True)),
            "mesh_poisson_depth"   : int(getattr(c, "mesh_level1_poisson_depth", 8)),
            "mesh_poisson_ptweight": float(getattr(c, "mesh_level1_poisson_ptweight", 0.5)),
        }

    def run(self, ctx: StageContext) -> None:
        import json
        from pathlib import Path

        import numpy as np
        import pyvista as pv
        from scipy import ndimage
        from scipy.spatial import cKDTree
        from sklearn.cluster import DBSCAN

        c = ctx.config
        stage_dir = Path(ctx.store.stage_dir(self.key))
        stage_dir.mkdir(parents=True, exist_ok=True)

        refined_world = np.asarray(ctx.require("refined_cluster"), dtype=np.float32)
        
        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        
        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = str(ctx.require("alpha_shell_path"))
        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        voxel = float(c.refine_voxel_size_A)
        pad = float(c.refine_roi_pad_A)
        atom_r = float(c.refine_atom_radius_A)
        keep_within_A = float(c.refine_keep_within_A)
        occ_close_iters = int(c.refine_occ_close_iters)
        void_open_iters = int(c.refine_void_open_iters)
        forbid_roi_boundary = bool(c.refine_forbid_roi_boundary)

        eps_c = float(c.dbscan_level1_coarse_eps_A)
        ms_c = int(c.dbscan_level1_coarse_min_samples)
        eps_r = float(c.dbscan_level1_refine_eps_A)
        ms_r = int(c.dbscan_level1_refine_min_samples)

        dbscan_max_points = int(c.refine_dbscan_max_points)
        dbscan_seed = int(c.refine_dbscan_seed)

        print(f"[{self.key}] voxel={voxel}Å, ROI pad={pad}Å, atom_r={atom_r}Å")

        refined_c0 = transform_points_to_C0(refined_world, ptc, constr).astype(np.float32)
        lo = refined_c0.min(axis=0) - pad
        hi = refined_c0.max(axis=0) + pad

        roi_obj = {
            "roi_id": "bbox_pad_stage50",
            "frame": "C0",
            "pad_A": float(pad),
            "lo": [float(x) for x in lo.tolist()],
            "hi": [float(x) for x in hi.tolist()],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
            "source": {"stage": "50_clustering", "artifact": "refined_cluster"},
        }
        (stage_dir / "roi_bbox_c0.json").write_text(json.dumps(roi_obj, indent=2))

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr).astype(np.float32)
        lo_sel = lo - atom_r
        hi_sel = hi + atom_r
        m_atoms = np.all((region_c0 >= lo_sel[None, :]) & (region_c0 <= hi_sel[None, :]), axis=1)
        atoms_roi_c0 = region_c0[m_atoms]
        
        print(f"[{self.key}] selected {atoms_roi_c0.shape[0]:,} atoms near ROI (ALL atoms, prevents interference)")

        grid = _make_bbox_grid(lo, hi, voxel)
        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))
        z_max = z_min + float(c.cylinder_height_A)

        cyl = self._cylinder_mask_bbox_grid(
            grid,
            radius_A=float(c.cylinder_radius_A),
            zmin_A=z_min,
            zmax_A=z_max,
        )

        occupied = occupancy_via_edt(atoms_roi_c0, grid, atom_radius_A=atom_r)

        if occ_close_iters > 0:
            occupied = ndimage.binary_closing(occupied, iterations=occ_close_iters)

        occupied = occupied | (~cyl)
        np.save(stage_dir / "occupied_mask_level_1.npy", occupied.astype(np.uint8))

        empty_mask = (~occupied) & cyl
        if empty_mask.sum() == 0:
            raise ValueError(f"[{self.key}] no empty voxels in ROI")

        empty_idx = np.argwhere(empty_mask)
        empty_pts_c0 = _voxel_centers_from_indices(grid, empty_idx).astype(np.float32)

        shell_world = pv.read(alpha_shell_path).triangulate()
        shell_c0 = shell_world.copy(deep=True)
        shell_c0.points = transform_points_to_C0(
            np.asarray(shell_world.points, dtype=np.float32), ptc, constr
        )

        sel = pv.PolyData(empty_pts_c0).select_enclosed_points(shell_c0, check_surface=watertight)
        inside_flags = (np.asarray(sel["SelectedPoints"], dtype=np.int8) == 1)

        void_mask = np.zeros_like(empty_mask, dtype=np.bool_)
        inside_idx = empty_idx[inside_flags]
        if inside_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] no empty voxels inside shell")
        
        void_mask[inside_idx[:, 0], inside_idx[:, 1], inside_idx[:, 2]] = True

        if forbid_roi_boundary:
            void_mask[0, :, :] = False
            void_mask[-1, :, :] = False
            void_mask[:, 0, :] = False
            void_mask[:, -1, :] = False
            void_mask[:, :, 0] = False
            void_mask[:, :, -1] = False

        if keep_within_A > 0.0:
            coarse_ijk = _points_to_ijk(grid, refined_c0)
            m_valid = _valid_ijk(grid, coarse_ijk)
            coarse_ijk = coarse_ijk[m_valid]
            if coarse_ijk.shape[0] > 0:
                seed = np.zeros(grid.shape, dtype=np.bool_)
                seed[coarse_ijk[:, 0], coarse_ijk[:, 1], coarse_ijk[:, 2]] = True
                dist_vox = ndimage.distance_transform_edt(~seed)
                r_vox = float(keep_within_A) / float(voxel)
                void_mask = void_mask & (dist_vox <= r_vox)

        if void_open_iters > 0:
            st = ndimage.generate_binary_structure(3, 1)
            void_mask = ndimage.binary_opening(void_mask, structure=st, iterations=void_open_iters)

        np.save(stage_dir / "void_mask_level_1.npy", void_mask.astype(np.uint8))

        st_er = ndimage.generate_binary_structure(3, 1)
        er = ndimage.binary_erosion(void_mask, structure=st_er, iterations=1)
        boundary_mask = void_mask & (~er)
        boundary_idx = np.argwhere(boundary_mask)
        
        if boundary_idx.shape[0] == 0:
            raise ValueError(f"[{self.key}] boundary extraction produced 0 voxels")

        boundary_pts_c0 = _voxel_centers_from_indices(grid, boundary_idx).astype(np.float32)
        boundary_pts_w = transform_points_from_C0(boundary_pts_c0, ptc, constr).astype(np.float32)

        print(f"[{self.key}] boundary points: {boundary_pts_w.shape[0]:,}")

        boundary_pts_c0_cap, boundary_pts_w_cap, idx_cap = self._maybe_cap(
            boundary_pts_c0, boundary_pts_w, dbscan_max_points, dbscan_seed
        )

        tree = cKDTree(refined_c0)

        t0 = time.perf_counter()
        db_coarse = DBSCAN(eps=eps_c, min_samples=ms_c, metric="euclidean", n_jobs=-1)
        labels_coarse = db_coarse.fit_predict(boundary_pts_c0_cap).astype(np.int32)
        dt0 = time.perf_counter() - t0

        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        d_coarse = stage_dir / "coarse"
        self._save_dbscan_pass(
            d_coarse,
            boundary_pts_w_cap,
            labels_coarse,
            eps=eps_c,
            min_samples=ms_c,
        )

        stats_coarse = self._cluster_stats(boundary_pts_c0_cap, labels_coarse, tree)
        best_coarse = self._choose_best_label(stats_coarse)
        
        if best_coarse == -1:
            raise ValueError(
                f"[{self.key}] coarse DBSCAN produced no clusters. "
                f"Try eps={eps_c} → {eps_c*1.5:.1f} or min_samples={ms_c} → {ms_c//2}"
            )

        m_best = labels_coarse == best_coarse
        pts_refine_c0 = boundary_pts_c0_cap[m_best]
        pts_refine_w = boundary_pts_w_cap[m_best]

        t1 = time.perf_counter()
        db_refine = DBSCAN(eps=eps_r, min_samples=ms_r, metric="euclidean", n_jobs=-1)
        labels_refine = db_refine.fit_predict(pts_refine_c0).astype(np.int32)
        dt1 = time.perf_counter() - t1

        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        d_refine = stage_dir / "refine"
        self._save_dbscan_pass(
            d_refine,
            pts_refine_w,
            labels_refine,
            eps=eps_r,
            min_samples=ms_r,
        )

        stats_refine = self._cluster_stats(pts_refine_c0, labels_refine, tree)
        best_refine = self._choose_best_label(stats_refine)
        
        if best_refine == -1:
            raise ValueError(
                f"[{self.key}] refine DBSCAN produced no clusters. "
                f"Try eps={eps_r} → {eps_r*1.5:.1f} or min_samples={ms_r} → {ms_r//2}"
            )

        final_mask = labels_refine == best_refine
        final_surface_w = pts_refine_w[final_mask].astype(np.float32)
        
        if final_surface_w.shape[0] == 0:
            raise ValueError(f"[{self.key}] final cluster is empty")

        np.save(stage_dir / "refined_surface_points_level_1.npy", final_surface_w)

        print(f"[{self.key}] creating selected void mask for voxel fallback...")
        
        final_surface_c0 = transform_points_to_C0(final_surface_w, ptc, constr).astype(np.float32)
        final_ijk = _points_to_ijk(grid, final_surface_c0)
        m_valid = _valid_ijk(grid, final_ijk)
        final_ijk = final_ijk[m_valid]
        
        boundary_selected = np.zeros(grid.shape, dtype=np.bool_)
        if final_ijk.shape[0] > 0:
            boundary_selected[final_ijk[:, 0], final_ijk[:, 1], final_ijk[:, 2]] = True
        
        from scipy.ndimage import binary_fill_holes, binary_dilation
        
        boundary_dilated = binary_dilation(boundary_selected, iterations=2)
        selected_mask = binary_fill_holes(boundary_dilated)
        selected_mask = selected_mask & void_mask
        
        selected_mask_path = stage_dir / "selected_void_component_mask_level_1.npy"
        np.save(selected_mask_path, selected_mask.astype(np.uint8))
        
        print(f"[{self.key}]   selected mask: {selected_mask.sum():,} voxels")
        
        ctx.inputs["selected_void_component_mask_level_1_path"] = str(selected_mask_path)
        ctx.inputs["grid_spec_level_1_path"] = str(stage_dir / "grid_spec_level_1.json")

        spec_obj = {
            "frame": "C0",
            "origin": [float(x) for x in grid.origin.tolist()],
            "voxel_size_A": float(grid.voxel_size),
            "shape": [int(x) for x in grid.shape],
            "transform": {
                "ptc": [float(x) for x in ptc.tolist()],
                "constriction": [float(x) for x in constr.tolist()]
            },
        }
        (stage_dir / "grid_spec_level_1.json").write_text(json.dumps(spec_obj, indent=2))

        diag = {
            "voxel_size_A": voxel,
            "roi_pad_A": pad,
            "atom_radius_A": atom_r,
            "keep_within_A": keep_within_A,
            "boundary_points_total": int(boundary_pts_c0.shape[0]),
            "boundary_points_used_for_dbscan": int(boundary_pts_c0_cap.shape[0]),
            "dbscan_coarse": {
                "eps_A": eps_c,
                "min_samples": ms_c,
                "best_label": int(best_coarse),
                "n_clusters": int(len(stats_coarse)),
                "clusters": stats_coarse[:25],
            },
            "dbscan_refine": {
                "eps_A": eps_r,
                "min_samples": ms_r,
                "best_label": int(best_refine),
                "n_clusters": int(len(stats_refine)),
                "clusters": stats_refine[:25],
            },
            "final_surface_points": int(final_surface_w.shape[0]),
        }
        (stage_dir / "dbscan_diagnostics.json").write_text(json.dumps(diag, indent=2))

        print(
            f"[{self.key}] DBSCAN refined-grid: "
            f"coarse={best_coarse} → refine={best_refine}, final_pts={final_surface_w.shape[0]:,}"
        )

        ctx.inputs["refined_cluster_surface"] = True
        ctx.inputs["refined_cluster"] = final_surface_w
        ctx.artifacts["refined_surface_points_level_1"] = str(stage_dir / "refined_surface_points_level_1.npy")

        if c.mesh_level1_enable:
            self._generate_mesh(ctx, final_surface_w, level_name="level_1")

    def _cylinder_mask_bbox_grid(
        self, grid: GridSpec, radius_A: float, zmin_A: float, zmax_A: float
    ) -> np.ndarray:
        nx, ny, nz = grid.shape
        v = float(grid.voxel_size)
        ox, oy, oz = grid.origin

        x = ox + np.arange(nx, dtype=np.float32) * v
        y = oy + np.arange(ny, dtype=np.float32) * v
        z = oz + np.arange(nz, dtype=np.float32) * v

        X, Y = np.meshgrid(x, y, indexing="ij")
        inside_r = (X * X + Y * Y) <= (radius_A * radius_A)
        inside_z = (z >= zmin_A) & (z <= zmax_A)
        return inside_r[:, :, None] & inside_z[None, None, :]

    def _maybe_cap(
        self, 
        points_c0: np.ndarray, 
        points_w: np.ndarray, 
        cap: int, 
        seed: int
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        n = points_c0.shape[0]
        if cap and n > cap:
            rng = np.random.default_rng(seed)
            idx = rng.choice(n, size=int(cap), replace=False)
            idx.sort()
            return points_c0[idx], points_w[idx], idx
        idx = np.arange(n, dtype=np.int64)
        return points_c0, points_w, idx

    def _cluster_stats(
        self, points_c0: np.ndarray, labels: np.ndarray, tree: cKDTree
    ) -> list[dict]:
        stats = []
        for lab in np.unique(labels):
            if lab == -1:
                continue
            m = labels == lab
            pts = points_c0[m]
            if pts.shape[0] == 0:
                continue
            d, _ = tree.query(pts, k=1)
            stats.append({
                "label": int(lab),
                "size": int(pts.shape[0]),
                "median_dist_to_stage50_A": float(np.median(d)),
                "p05_dist_to_stage50_A": float(np.percentile(d, 5)),
                "p95_dist_to_stage50_A": float(np.percentile(d, 95)),
            })
        stats.sort(key=lambda x: (x["median_dist_to_stage50_A"], -x["size"]))
        return stats

    def _choose_best_label(self, stats: list[dict]) -> int:
        return int(stats[0]["label"]) if stats else -1

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        eps: float,
        min_samples: int,
    ) -> None:
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        clusters_from_labels_func = __import__('ribctl.lib.npet2.backends.clustering_io', fromlist=['clusters_from_labels']).clusters_from_labels
        clusters = clusters_from_labels_func(pts, labels)
        for cid, cpts in clusters.items():
            if cid == -1:
                continue
            if cpts.shape[0] > 0:
                np.save(pass_dir / f"cluster_id{cid}.npy", cpts.astype(np.float32))

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        import time
        import json
        from ribctl.lib.npet.kdtree_approach import transform_points_from_C0
        from ribctl.lib.npet2.backends.meshing import (
            mesh_from_binary_volume,
            clip_mesh_to_atom_clearance,
            save_mesh_with_ascii,
        )

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        print(f"[{self.key}] generating mesh for {level_name}...")

        mask_path = stage_dir / "selected_void_component_mask_level_1.npy"
        if not mask_path.exists():
            mask_path = stage_dir / "void_mask_level_1.npy"
        if not mask_path.exists():
            print(f"[{self.key}] no void mask found for {level_name}, skipping mesh")
            return

        mask = np.load(mask_path).astype(bool)
        spec = json.loads((stage_dir / "grid_spec_level_1.json").read_text())
        origin = np.asarray(spec["origin"], dtype=np.float32)
        voxel = float(spec["voxel_size_A"])
        ptc = np.asarray(spec["transform"]["ptc"], dtype=np.float32)
        constr = np.asarray(spec["transform"]["constriction"], dtype=np.float32)

        t0 = time.perf_counter()
        try:
            surf_c0, pre_smooth_c0 = mesh_from_binary_volume(
                mask, origin, voxel,
                gaussian_sigma_voxels=c.mesh_level1_gaussian_sigma,
                smooth_method=c.mesh_smooth_method,
                smooth_iters=c.mesh_level1_smooth_iters,
                taubin_pass_band=c.mesh_taubin_pass_band,
                fill_holes_size=c.mesh_fill_holes_A,
            )
        except ValueError as e:
            print(f"[{self.key}] MC mesh failed for {level_name}: {e}")
            return

        def _to_world(mesh_c0: pv.PolyData) -> pv.PolyData:
            pts_w = transform_points_from_C0(
                np.asarray(mesh_c0.points, dtype=np.float32), ptc, constr
            ).astype(np.float32)
            m = mesh_c0.copy(deep=True)
            m.points = pts_w
            return m

        pre_smooth_w = _to_world(pre_smooth_c0)
        pre_smooth_path = stage_dir / f"mesh_{level_name}_pre_smooth.ply"
        save_mesh_with_ascii(pre_smooth_w, pre_smooth_path, tag=f"{level_name}-pre-smooth")

        surf_w = _to_world(surf_c0)

        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)
        surf_w = clip_mesh_to_atom_clearance(surf_w, region_xyz, min_clearance_A=c.mesh_atom_clearance_A)

        dt = time.perf_counter() - t0
        is_watertight = surf_w.is_manifold and surf_w.n_open_edges == 0
        print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
              f"{surf_w.n_faces:,} faces, watertight={is_watertight}")

        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        save_mesh_with_ascii(surf_w, mesh_path, tag=level_name)

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name, "method": "marching_cubes_taubin",
                  "watertight": is_watertight, "voxel_size_A": voxel},
        )
        print(f"[{self.key}] mesh saved: {mesh_path}")

        ctx.inputs["level_1_mesh_path"] = str(mesh_path)
        ctx.inputs["level_1_mesh_watertight"] = is_watertight

    def _clip_mesh_to_atoms(
        self,
        mesh: pv.PolyData,
        atom_xyz: np.ndarray,
        min_clearance_A: float = 1.5,
    ) -> pv.PolyData:
        """
        Push any mesh vertices that are closer than min_clearance_A to any atom
        back along the atom->vertex direction until they are at min_clearance_A.
        
        This prevents Poisson/smoothing overshoot from eating into atom-occupied space.
        """
        from scipy.spatial import cKDTree

        tree = cKDTree(atom_xyz)
        pts = np.asarray(mesh.points, dtype=np.float64)

        dist, idx = tree.query(pts, k=1)
        violating = dist < min_clearance_A
        n_violations = int(violating.sum())

        if n_violations == 0:
            print(f"[{self.key}]   atom clearance check: all vertices OK (min_clearance={min_clearance_A}A)")
            return mesh

        nearest_atom = atom_xyz[idx[violating]]
        direction = pts[violating] - nearest_atom
        norms = np.linalg.norm(direction, axis=1, keepdims=True)
        norms = np.maximum(norms, 1e-8)
        direction = direction / norms

        pts[violating] = nearest_atom + direction * min_clearance_A

        result = mesh.copy()
        result.points = pts.astype(np.float32)

        print(f"[{self.key}]   atom clearance check: pushed {n_violations:,} vertices "
            f"(of {pts.shape[0]:,}) to min {min_clearance_A}A clearance")

        return result
```

ribctl/lib/npet2/stages/legacy_minimal.py
```py
# ribctl/lib/npet2/stages/legacy_minimal.py

from __future__ import annotations
import json
from pathlib import Path
import time
from typing import Any, Dict, List, Tuple
import numpy as np
import pyvista as pv
import open3d as o3d

from ribctl.lib.npet2.backends.grid_occupancy import (
    connected_components_3d,
    occupancy_via_edt,
)
from ribctl.lib.npet2.backends.meshing import save_mesh_with_ascii
from ribctl.lib.npet2.core.cache import StageCacheKey
from ribctl.lib.npet2.core.pipeline import Stage
from ribctl.lib.npet2.core.ribosome_types import RibosomeProfile
from ribctl.lib.npet2.core.structure_selection import (
    intersect_with_first_assembly,
    ribosome_wall_auth_asym_ids,
    tunnel_debris_chains,
    atom_inclusion_policy,
)
from ribctl.lib.npet2.core.types import StageContext, ArtifactType

from scipy import ndimage

from ribctl.lib.npet.alphalib import (
    cif_to_point_cloud,
    fast_normal_estimation,
    quick_surface_points,
    validate_mesh_pyvista,
)
from ribctl.lib.npet.kdtree_approach import (
    apply_poisson_reconstruction,
    filter_residues_parallel,
    transform_points_to_C0,
    transform_points_from_C0,
    create_point_cloud_mask,
    DBSCAN_capture,
    DBSCAN_pick_largest_cluster,
    estimate_normals,
)
from ribctl.lib.npet2.stages.grid_refine import (
    _make_bbox_grid,
    _points_to_ijk,
    _valid_ijk,
    _voxel_centers_from_indices,
)


def _residues_from_chain_ids(structure, chain_ids: set[str]):
    model = structure[0]
    residues = []
    for cid in chain_ids:
        if cid not in model:
            continue
        chain = model[cid]
        for r in chain.get_residues():
            if len(getattr(r, "child_list", [])) == 0:
                continue
            residues.append(r)
    return residues


def _pick_tunnel_cluster(
    clusters: dict[int, list],
    constr: np.ndarray,
) -> tuple[np.ndarray, int]:
    """Pick the cluster whose points are closest to the constriction site."""
    constr = np.asarray(constr, dtype=np.float32).reshape(1, 3)
    best_id = -1
    best_dist = float("inf")
    for cid, pts_list in clusters.items():
        if cid == -1:
            continue
        pts = np.asarray(pts_list, dtype=np.float32)
        if pts.shape[0] == 0:
            continue
        dists = np.linalg.norm(pts - constr, axis=1)
        min_dist = float(dists.min())
        if min_dist < best_dist:
            best_dist = min_dist
            best_id = cid
    if best_id == -1:
        raise ValueError("No valid clusters found")
    print(f"  [cluster_select] picked cluster {best_id} "
          f"(n={len(clusters[best_id]):,}, dist_to_constriction={best_dist:.1f}A)")
    return np.asarray(clusters[best_id], dtype=np.float32), best_id

def _get_biopython_structure(ctx: StageContext):
    """Get biopython structure from ctx, or parse mmcif on demand."""
    bs = ctx.inputs.get("biopython_structure")
    if bs is not None:
        return bs
    from Bio.PDB.MMCIFParser import FastMMCIFParser
    mmcif_path = ctx.require("mmcif_path")
    bs = FastMMCIFParser(QUIET=True).get_structure(ctx.rcsb_id, mmcif_path)
    ctx.inputs["biopython_structure"] = bs
    return bs


class Stage20ExteriorShell(Stage):
    key = "20_exterior_shell"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "d3d_alpha": c.alpha_d3d_alpha,
            "d3d_tol": c.alpha_d3d_tol,
            "d3d_offset": c.alpha_d3d_offset,
            "kdtree_radius": c.alpha_kdtree_radius,
            "max_nn": c.alpha_max_nn,
            "tangent_k": c.alpha_tangent_planes_k,
            "poisson_depth": c.alpha_poisson_depth,
            "poisson_ptweight": c.alpha_poisson_ptweight,
            "fill_holes": c.alpha_fill_holes,
        }

    def run(self, ctx: StageContext) -> None:
        stage_cache = ctx.require("stage_cache")
        inputs_fp = ctx.require("inputs_fp")
        params = self.params(ctx)

        key = StageCacheKey(
            stage=self.key,
            inputs_fp={"structure": inputs_fp["structure"]},
            params=params,
            impl_version="v1",
        )

        stage_dir = ctx.store.stage_dir(self.key)
        cached_files = [
            "alpha_shell.ply",
            "alpha_shell_quality.json",
            "alpha_normals.ply",
            "alpha_surface_points.npy",
            "ribosome_ptcloud.npy",
        ]

        if stage_cache.has(
            key, required=["alpha_shell.ply", "alpha_shell_quality.json"]
        ):
            stage_cache.copy_into(key, stage_dir, cached_files)

            quality = json.loads((stage_dir / "alpha_shell_quality.json").read_text())
            ctx.inputs["alpha_shell_path"] = str(stage_dir / "alpha_shell.ply")
            ctx.inputs["alpha_shell_watertight"] = bool(
                quality.get("watertight", False)
            )

            ctx.store.register_file(
                name="alpha_shell_mesh",
                stage=self.key,
                type=ArtifactType.PLY_MESH,
                path=stage_dir / "alpha_shell.ply",
            )
            ctx.store.register_file(
                name="alpha_shell_quality",
                stage=self.key,
                type=ArtifactType.JSON,
                path=stage_dir / "alpha_shell_quality.json",
            )
            return

        c = ctx.config
        profile: RibosomeProfile = ctx.require("profile")
        cifpath = Path(ctx.require("mmcif_path"))

        ptcloud_path = stage_dir / "ribosome_ptcloud.npy"
        surface_pts_path = stage_dir / "alpha_surface_points.npy"
        normals_pcd_path = stage_dir / "alpha_normals.ply"
        mesh_path = stage_dir / "alpha_shell.ply"
        quality_path = stage_dir / "alpha_shell_quality.json"

        wall = ribosome_wall_auth_asym_ids(
            profile,
            exclude_trna=bool(getattr(ctx.config, "occupancy_exclude_trna", True)),
            extra_exclude=tunnel_debris_chains(ctx.rcsb_id, profile),
        )
        wall = intersect_with_first_assembly(profile, wall)

        ptcloud = cif_to_point_cloud(str(cifpath), sorted(wall), do_atoms=True)

        np.save(ptcloud_path, ptcloud)
        ctx.store.register_file(
            name="ribosome_ptcloud",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=ptcloud_path,
        )

        surface_pts = quick_surface_points(
            ptcloud, c.alpha_d3d_alpha, c.alpha_d3d_tol, c.alpha_d3d_offset
        ).astype(np.float32)
        np.save(surface_pts_path, surface_pts)
        ctx.store.register_file(
            name="alpha_surface_points",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=surface_pts_path,
        )

        normal_estimated_pcd = fast_normal_estimation(
            surface_pts, c.alpha_kdtree_radius, c.alpha_max_nn, c.alpha_tangent_planes_k
        )

        center = normal_estimated_pcd.get_center()
        normal_estimated_pcd.orient_normals_towards_camera_location(
            camera_location=center
        )
        normal_estimated_pcd.normals = o3d.utility.Vector3dVector(
            -np.asarray(normal_estimated_pcd.normals)
        )

        o3d.io.write_point_cloud(str(normals_pcd_path), normal_estimated_pcd)
        ctx.store.register_file(
            name="alpha_normals_pcd",
            stage=self.key,
            type=ArtifactType.PLY_PCD,
            path=normals_pcd_path,
        )

        apply_poisson_reconstruction(
            str(normals_pcd_path),
            mesh_path,
            recon_depth=c.alpha_poisson_depth,
            recon_pt_weight=c.alpha_poisson_ptweight,
        )

        mesh = pv.read(mesh_path)
        mesh = mesh.fill_holes(c.alpha_fill_holes)
        mesh = mesh.connectivity(largest=True).triangulate()
        mesh.save(mesh_path)

        watertight = validate_mesh_pyvista(mesh)

        quality = {
            "watertight": bool(watertight),
            "n_points": int(mesh.n_points),
            "n_faces": int(mesh.n_faces),
            "open_edges": int(mesh.n_open_edges),
            "is_manifold": bool(mesh.is_manifold),
            "bounds": list(mesh.bounds),
        }
        quality_path.write_text(json.dumps(quality, indent=2))
        ctx.store.register_file(
            name="alpha_shell_quality",
            stage=self.key,
            type=ArtifactType.JSON,
            path=quality_path,
        )

        ctx.store.register_file(
            name="alpha_shell_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
        )
        ctx.inputs["alpha_shell_path"] = str(mesh_path)
        ctx.inputs["alpha_shell_watertight"] = bool(watertight)
        if watertight:
            stage_cache.put_from(key, stage_dir, cached_files)

class Stage30RegionAtoms(Stage):
    key = "30_region_atoms"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {"radius_A": c.cylinder_radius_A, "height_A": c.cylinder_height_A}

    def run(self, ctx: StageContext) -> None:
        c = ctx.config
        profile: RibosomeProfile = ctx.require("profile")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        policy = atom_inclusion_policy(profile, c, ctx.rcsb_id)

        occ_chain_ids = policy["wall_chain_ids"]
        seed_chain_ids = set(occ_chain_ids)

        structure = _get_biopython_structure(ctx)

        residues_seed = _residues_from_chain_ids(structure, seed_chain_ids)
        residues_occ = _residues_from_chain_ids(structure, occ_chain_ids)

        residues_seed = filter_residues_parallel(
            residues=residues_seed,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
        )
        residues_occ = filter_residues_parallel(
            residues=residues_occ,
            base_point=ptc,
            axis_point=constr,
            radius=c.cylinder_radius_A,
            height=c.cylinder_height_A,
            max_workers=1,
            chunk_size=5000,
            z_min=z_min,
        )

        seed_points = np.asarray(
            [atom.get_coord() for r in residues_seed for atom in r.child_list],
            dtype=np.float32,
        )
        occ_points = np.asarray(
            [atom.get_coord() for r in residues_occ for atom in r.child_list],
            dtype=np.float32,
        )

        stage_dir = ctx.store.stage_dir(self.key)

        out_seed = stage_dir / "region_atom_xyz.npy"
        np.save(out_seed, seed_points)
        ctx.store.register_file(
            name="region_atom_xyz",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_seed,
            meta={
                "n": int(seed_points.shape[0]),
                "note": "seed atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz"] = seed_points

        out_occ = stage_dir / "region_atom_xyz_occ.npy"
        np.save(out_occ, occ_points)
        ctx.store.register_file(
            name="region_atom_xyz_occ",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=out_occ,
            meta={
                "n": int(occ_points.shape[0]),
                "note": "occupancy atoms (walls-only chains)",
            },
        )
        ctx.inputs["region_atom_xyz_occ"] = occ_points

        policy_record = {
            "occupancy_chain_mode": getattr(c, "occupancy_chain_mode", "walls_only"),
            "exclude_trna": bool(getattr(c, "occupancy_exclude_trna", True)),
            "wall_chain_ids": sorted(occ_chain_ids),
            "excluded_chains": {k: v for k, v in policy["reasons"].items()},
            "policy_summary": (
                "INCLUDED: ribosomal proteins + rRNAs (with modified residues). "
                "EXCLUDED: waters, ions, nonpolymer ligands, tRNAs, debris chains."
            ),
        }
        (stage_dir / "atom_selection_policy.json").write_text(
            json.dumps(policy_record, indent=2)
        )

        print(
            f"[{self.key}] seed_atoms={seed_points.shape[0]:,} "
            f"occ_atoms={occ_points.shape[0]:,} occ_chains={len(occ_chain_ids)}"
        )
        if policy["reasons"]:
            excluded_summary = ", ".join(
                f"{k}({v})" for k, v in sorted(policy["reasons"].items())
            )
            print(f"[{self.key}] excluded: {excluded_summary}")

class Stage40EmptySpace(Stage):
    key = "40_empty_space"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "grid_levels": [
                {
                    "name": gl.name,
                    "voxel_size_A": gl.voxel_size_A,
                    "backend": gl.occupancy_backend,
                    "atom_radius_mode": getattr(gl, "atom_radius_mode", "uniform"),
                    "uniform_atom_radius_A": getattr(gl, "uniform_atom_radius_A", None),
                }
                for gl in c.grid_levels
            ],
            "radius_A": c.cylinder_radius_A,
            "height_A": c.cylinder_height_A,
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import numpy as np
        import pyvista as pv

        from ribctl.lib.npet.kdtree_approach import (
            create_point_cloud_mask,
            transform_points_from_C0,
            transform_points_to_C0,
        )

        from ribctl.lib.npet2.backends.grid_occupancy import (
            make_cylinder_grid,
            cylinder_mask,
            occupancy_via_edt,
            empty_points_from_mask,
            save_grid_npy,
            get_occupied_voxel_centers,
        )

        c = ctx.config

        z_min = float(ctx.inputs.get("cylinder_z_min", 0.0))

        # Use occ atoms for occupancy (prevents mesh interference)
        # region_xyz = np.asarray(ctx.require("region_atom_xyz_all"), dtype=np.float32)
        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)

        # Use filtered atoms for clustering seed reference
        region_xyz_filtered = np.asarray(
            ctx.require("region_atom_xyz"), dtype=np.float32
        )

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)
        alpha_shell_path = ctx.require("alpha_shell_path")

        region_c0 = transform_points_to_C0(region_xyz, ptc, constr)

        shell = pv.read(alpha_shell_path)
        if not isinstance(shell, pv.PolyData):
            shell = shell.extract_surface()
        if not shell.is_all_triangles:
            shell = shell.triangulate()

        watertight = bool(ctx.inputs.get("alpha_shell_watertight", True))

        stage_dir = ctx.store.stage_dir(self.key)

        clip_note = {
            "alpha_shell_path": str(alpha_shell_path),
            "alpha_shell_watertight": watertight,
            "clipping_mode": "select_enclosed_points(check_surface=True)"
            if watertight
            else "select_enclosed_points(check_surface=False) [fallback; shell not watertight]",
            "occupancy_atoms": "region_atom_xyz_all (includes ALL chains)",
        }
        p_clip_note = stage_dir / "clipping_note.json"
        p_clip_note.write_text(json.dumps(clip_note, indent=2))
        ctx.store.register_file(
            name="clipping_note",
            stage=self.key,
            type=ArtifactType.JSON,
            path=p_clip_note,
        )

        last_empty = None

        for gl in c.grid_levels:
            backend = gl.occupancy_backend

            if backend == "legacy_kdtree":
                mask, (x, y, z) = create_point_cloud_mask(
                    region_c0,
                    radius=c.cylinder_radius_A,
                    height=c.cylinder_height_A,
                    voxel_size=gl.voxel_size_A,
                    radius_around_point=gl.uniform_atom_radius_A,
                    z_min=z_min,
                )

                idx = np.where(~mask)
                empty_c0 = np.column_stack((x[idx[0]], y[idx[1]], z[idx[2]])).astype(
                    np.float32
                )

            elif backend == "edt":
                grid = make_cylinder_grid(
                    radius_A=float(c.cylinder_radius_A),
                    height_A=float(c.cylinder_height_A),
                    voxel_A=float(gl.voxel_size_A),
                    z_min=z_min,
                )

                occ = occupancy_via_edt(
                    region_c0,
                    grid,
                    atom_radius_A=float(gl.uniform_atom_radius_A),
                )

                cyl2d = cylinder_mask(grid, radius_A=float(c.cylinder_radius_A))
                cyl = np.broadcast_to(cyl2d, grid.shape)

                occ = occ | (~cyl)

                empty_mask = ~occ
                empty_c0 = empty_points_from_mask(grid, empty_mask & cyl)

                save_grid_npy(
                    grid, occ, stage_dir / f"occupancy_grid_{gl.name}", compress=False
                )
                save_grid_npy(
                    grid,
                    (empty_mask & cyl),
                    stage_dir / f"empty_mask_{gl.name}",
                    compress=False,
                )

                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"occupancy_grid_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"occupancy_grid_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"occupancy_grid_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_data",
                    stage=self.key,
                    type=ArtifactType.NUMPY,
                    path=stage_dir / f"empty_mask_{gl.name}_data.npy",
                    meta={"voxel_size_A": gl.voxel_size_A, "shape": list(grid.shape)},
                )
                ctx.store.register_file(
                    name=f"empty_mask_{gl.name}_spec",
                    stage=self.key,
                    type=ArtifactType.JSON,
                    path=stage_dir / f"empty_mask_{gl.name}_spec.json",
                    meta={"voxel_size_A": gl.voxel_size_A},
                )

                try:
                    occ_centers_c0 = get_occupied_voxel_centers(grid, occ).astype(
                        np.float32
                    )
                    occ_centers_world = transform_points_from_C0(
                        occ_centers_c0, ptc, constr
                    ).astype(np.float32)
                    p_occ = stage_dir / f"occupied_voxels_{gl.name}.npy"
                    np.save(p_occ, occ_centers_world)
                    ctx.store.register_file(
                        name=f"occupied_voxels_{gl.name}",
                        stage=self.key,
                        type=ArtifactType.NUMPY,
                        path=p_occ,
                        meta={
                            "voxel_size_A": gl.voxel_size_A,
                            "n": int(occ_centers_world.shape[0]),
                        },
                    )
                except Exception:
                    pass

            else:
                raise ValueError(
                    f"Grid level {gl.name}: unsupported backend {backend} "
                    f"(supported: legacy_kdtree, edt)"
                )

            empty_world = transform_points_from_C0(empty_c0, ptc, constr).astype(
                np.float32
            )

            p_pre = stage_dir / f"empty_points_{gl.name}_preclip.npy"
            np.save(p_pre, empty_world)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}_preclip",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=p_pre,
                meta={"voxel_size_A": gl.voxel_size_A, "n": int(empty_world.shape[0])},
            )

            if empty_world.shape[0] == 0:
                inside = empty_world
            else:
                pts_poly = pv.PolyData(empty_world)
                sel = pts_poly.select_enclosed_points(shell, check_surface=watertight)
                inside = empty_world[sel["SelectedPoints"] == 1].astype(np.float32)

            out = stage_dir / f"empty_points_{gl.name}.npy"
            np.save(out, inside)
            ctx.store.register_file(
                name=f"empty_points_{gl.name}",
                stage=self.key,
                type=ArtifactType.NUMPY,
                path=out,
                meta={
                    "voxel_size_A": gl.voxel_size_A,
                    "n": int(inside.shape[0]),
                    "backend": backend,
                    "alpha_shell_watertight": watertight,
                },
            )

            ctx.inputs[f"empty_points_{gl.name}"] = inside
            last_empty = inside

        ctx.inputs["empty_points"] = last_empty

class Stage50Clustering(Stage):
    """
    DBSCAN clustering on level_0 (coarse grid, typically 1.0A).
    
    Two-pass strategy:
      1. Coarse DBSCAN: merge regions, bridge gaps
      2. Refine DBSCAN: tighten on largest cluster from pass 1
    
    Cluster selection uses axial proximity to the PTC-Constriction axis
    rather than raw point count, which prevents the inter-subunit space
    from being picked over the actual tunnel.
    
    Optionally generates a mesh from the refined cluster.
    
    Outputs:
      stage/50_clustering/
        coarse/{points.npy, labels.npy, cluster_*.npy, index.json}
        refine/{points.npy, labels.npy, cluster_*.npy, index.json}
        largest_cluster.npy
        refined_cluster.npy
        mesh_level_0.ply (if enabled)
    """
    key = "50_clustering"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        c = ctx.config
        return {
            "coarse_eps_A": c.dbscan_level0_coarse_eps_A,
            "coarse_min_samples": c.dbscan_level0_coarse_min_samples,
            "refine_eps_A": c.dbscan_level0_refine_eps_A,
            "refine_min_samples": c.dbscan_level0_refine_min_samples,
            "mesh_enable": bool(getattr(c, "mesh_level0_enable", True)),
            "cluster_selection": "axial_proximity",
        }

    def run(self, ctx: StageContext) -> None:
        import json
        import time

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)

        # Input: empty points from Stage40
        empty_pts = np.asarray(ctx.require("empty_points"), dtype=np.float32)
        if empty_pts.ndim != 2 or empty_pts.shape[1] != 3 or empty_pts.shape[0] == 0:
            raise ValueError(f"[{self.key}] empty_points must be (N,3) and non-empty, got {empty_pts.shape}")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        print(f"[{self.key}] empty_points n={empty_pts.shape[0]:,}")

        # -----------------------
        # PASS 1: Coarse DBSCAN
        # -----------------------
        t0 = time.perf_counter()
        db_coarse, clusters_coarse = DBSCAN_capture(
            empty_pts, 
            c.dbscan_level0_coarse_eps_A, 
            c.dbscan_level0_coarse_min_samples
        )
        labels_coarse = np.asarray(db_coarse.labels_, dtype=np.int32)
        dt0 = time.perf_counter() - t0
        
        n_clusters_coarse = len(set(labels_coarse.tolist())) - (1 if -1 in labels_coarse else 0)
        print(f"[{self.key}] coarse DBSCAN: {dt0:.2f}s, {n_clusters_coarse} clusters")

        self._save_dbscan_pass(
            stage_dir / "coarse",
            empty_pts,
            labels_coarse,
            clusters_coarse,
            eps=c.dbscan_level0_coarse_eps_A,
            min_samples=c.dbscan_level0_coarse_min_samples,
        )

        # Pick tunnel cluster from coarse (axial proximity, not largest)


        largest, largest_id = _pick_tunnel_cluster(clusters_coarse, constr)
        if largest.shape[0] == 0:
            raise ValueError(f"[{self.key}] tunnel cluster is empty")

        p_largest = stage_dir / "largest_cluster.npy"
        np.save(p_largest, largest)
        ctx.store.register_file(
            name="largest_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_largest,
            meta={"cluster_id": int(largest_id), "n": int(largest.shape[0]),
                  "selection": "axial_proximity"},
        )

        # -----------------------
        # PASS 2: Refine DBSCAN
        # -----------------------
        t1 = time.perf_counter()
        db_refine, clusters_refine = DBSCAN_capture(
            largest,
            c.dbscan_level0_refine_eps_A,
            c.dbscan_level0_refine_min_samples,
        )
        labels_refine = np.asarray(db_refine.labels_, dtype=np.int32)
        dt1 = time.perf_counter() - t1
        
        n_clusters_refine = len(set(labels_refine.tolist())) - (1 if -1 in labels_refine else 0)
        print(f"[{self.key}] refine DBSCAN: {dt1:.2f}s, {n_clusters_refine} clusters")

        self._save_dbscan_pass(
            stage_dir / "refine",
            largest,
            labels_refine,
            clusters_refine,
            eps=c.dbscan_level0_refine_eps_A,
            min_samples=c.dbscan_level0_refine_min_samples,
        )

        # Pick tunnel cluster from refine


        refined, refined_id = _pick_tunnel_cluster(clusters_refine, constr)
        if refined.shape[0] == 0:
            raise ValueError(f"[{self.key}] refined cluster is empty")

        p_refined = stage_dir / "refined_cluster.npy"
        np.save(p_refined, refined)
        ctx.store.register_file(
            name="refined_cluster",
            stage=self.key,
            type=ArtifactType.NUMPY,
            path=p_refined,
            meta={"cluster_id": int(refined_id), "n": int(refined.shape[0]),
                  "selection": "axial_proximity"},
        )

        # Output for downstream
        ctx.inputs["refined_cluster"] = refined
        ctx.inputs["largest_cluster"] = largest

        print(f"[{self.key}] winner: coarse={largest.shape[0]:,} -> refine={refined.shape[0]:,}")

        # -----------------------
        # Optional: Mesh level_0
        # -----------------------
        if c.mesh_level0_enable:
            self._generate_mesh(ctx, refined, level_name="level_0")

    def _save_dbscan_pass(
        self,
        pass_dir: Path,
        pts: np.ndarray,
        labels: np.ndarray,
        clusters_dict: dict[int, list],
        eps: float,
        min_samples: int,
    ) -> None:
        """Save DBSCAN pass artifacts."""
        pass_dir.mkdir(parents=True, exist_ok=True)

        np.save(pass_dir / "points.npy", pts.astype(np.float32))
        np.save(pass_dir / "labels.npy", labels.astype(np.int32))

        counts = {}
        for lab in np.unique(labels):
            counts[int(lab)] = int((labels == lab).sum())

        index = {
            "pass": pass_dir.name,
            "eps_A": float(eps),
            "min_samples": int(min_samples),
            "n_points": int(pts.shape[0]),
            "labels": counts,
        }
        (pass_dir / "index.json").write_text(json.dumps(index, indent=2))

        for lab, plist in clusters_dict.items():
            lab = int(lab)
            if lab == -1:
                continue
            arr = np.asarray(plist, dtype=np.float32)
            if arr.size > 0:
                np.save(pass_dir / f"cluster_id{lab}.npy", arr)

    def _generate_mesh(self, ctx: StageContext, points: np.ndarray, level_name: str) -> None:
        import time
        from ribctl.lib.npet.kdtree_approach import transform_points_to_C0, transform_points_from_C0
        from ribctl.lib.npet2.backends.meshing import (
            mesh_from_binary_volume,
            voxelize_points,
            clip_mesh_to_atom_clearance,
            save_mesh_with_ascii,
        )

        c = ctx.config
        stage_dir = ctx.store.stage_dir(self.key)
        print(f"[{self.key}] generating mesh for {level_name}...")

        ptc = np.asarray(ctx.require("ptc_xyz"), dtype=np.float32)
        constr = np.asarray(ctx.require("constriction_xyz"), dtype=np.float32)

        pts_c0 = transform_points_to_C0(points, ptc, constr).astype(np.float32)
        voxel = 1.0

        t0 = time.perf_counter()
        mask, origin = voxelize_points(pts_c0, voxel_size=voxel, pad_voxels=2)

        try:
            surf_c0, pre_smooth_c0 = mesh_from_binary_volume(
                mask, origin, voxel,
                gaussian_sigma_voxels=c.mesh_level0_gaussian_sigma,
                smooth_method=c.mesh_smooth_method,
                smooth_iters=c.mesh_level0_smooth_iters,
                taubin_pass_band=c.mesh_taubin_pass_band,
                fill_holes_size=c.mesh_fill_holes_A,
            )
        except ValueError as e:
            print(f"[{self.key}] MC mesh failed for {level_name}: {e}")
            return

        # Transform both meshes to world coordinates
        def _to_world(mesh_c0: pv.PolyData) -> pv.PolyData:
            pts_w = transform_points_from_C0(
                np.asarray(mesh_c0.points, dtype=np.float32), ptc, constr
            ).astype(np.float32)
            m = mesh_c0.copy(deep=True)
            m.points = pts_w
            return m

        pre_smooth_w = _to_world(pre_smooth_c0)
        pre_smooth_path = stage_dir / f"mesh_{level_name}_pre_smooth.ply"
        save_mesh_with_ascii(pre_smooth_w, pre_smooth_path, tag=f"{level_name}-pre-smooth")

        surf_w = _to_world(surf_c0)

        region_xyz = np.asarray(ctx.require("region_atom_xyz_occ"), dtype=np.float32)
        surf_w = clip_mesh_to_atom_clearance(surf_w, region_xyz, min_clearance_A=c.mesh_atom_clearance_A)

        dt = time.perf_counter() - t0
        print(f"[{self.key}]   MC mesh: {dt:.2f}s, {surf_w.n_points:,} pts, "
              f"{surf_w.n_faces:,} faces, watertight={surf_w.is_manifold and surf_w.n_open_edges == 0}")

        mesh_path = stage_dir / f"mesh_{level_name}.ply"
        save_mesh_with_ascii(surf_w, mesh_path, tag=level_name)

        ctx.store.register_file(
            name=f"mesh_{level_name}",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"level": level_name, "method": "marching_cubes_taubin"},
        )
        print(f"[{self.key}] mesh saved: {mesh_path}")

class Stage70MeshValidate(Stage):
    key = "70_mesh_validate"

    def params(self, ctx: StageContext) -> Dict[str, Any]:
        return {}

    def run(self, ctx: StageContext) -> None:
        import json
        import shutil

        stage_dir = ctx.store.stage_dir(self.key)
        mesh_path = stage_dir / "npet2_tunnel_mesh.ply"

        def _mesh_stats(m: pv.PolyData) -> dict:
            return {
                "n_points": int(m.n_points),
                "n_faces": int(m.n_faces),
                "open_edges": int(m.n_open_edges),
                "is_manifold": bool(m.is_manifold),
                "bounds": [float(x) for x in m.bounds],
            }

        # Try level_1 first (higher detail), fall back to level_0
        chosen_src = None
        chosen_label = None

        l1_path = ctx.inputs.get("level_1_mesh_path")
        if l1_path and Path(l1_path).exists():
            m = pv.read(l1_path)
            if m.is_manifold and m.n_open_edges == 0 and m.n_points > 0:
                chosen_src = l1_path
                chosen_label = "level_1"
                print(f"[{self.key}] using level_1 mesh (0.5A grid)")

        if chosen_src is None:
            # Look for level_0
            l0_path = ctx.store.run_dir / "stage" / "50_clustering" / "mesh_level_0.ply"
            if l0_path.exists():
                m = pv.read(str(l0_path))
                if m.n_points > 0:
                    chosen_src = str(l0_path)
                    chosen_label = "level_0"
                    print(f"[{self.key}] falling back to level_0 mesh (1.0A grid)")

        if chosen_src is None:
            raise ValueError(f"[{self.key}] no valid mesh found from any stage")

        shutil.copy2(chosen_src, mesh_path)
        final = pv.read(str(mesh_path))

        # Save final mesh as both binary and ASCII
        save_mesh_with_ascii(final, mesh_path, tag="final")

        st = _mesh_stats(final)
        watertight = final.is_manifold and final.n_open_edges == 0
        print(f"[{self.key}] final mesh: {st}, watertight={watertight}")

        if not watertight:
            raise ValueError(f"[{self.key}] final mesh is not watertight")

        ctx.store.register_file(
            name="tunnel_mesh",
            stage=self.key,
            type=ArtifactType.PLY_MESH,
            path=mesh_path,
            meta={"watertight": True, "source": chosen_label},
        )
        ctx.inputs["tunnel_mesh_path"] = str(mesh_path)

        # Also copy final meshes to run root for convenience
        root_mesh = ctx.store.run_dir / "tunnel_mesh.ply"
        root_mesh_ascii = ctx.store.run_dir / "tunnel_mesh_ascii.ply"
        shutil.copy2(str(mesh_path), str(root_mesh))
        ascii_src = mesh_path.parent / f"{mesh_path.stem}_ascii.ply"
        if ascii_src.exists():
            shutil.copy2(str(ascii_src), str(root_mesh_ascii))

        self._copy_comparison_meshes(ctx, stage_dir)

    def _copy_comparison_meshes(self, ctx, stage_dir):
        import shutil

        run_root = ctx.store.run_dir

        for stage_name, level, voxel in [
            ("50_clustering", "level_0", 1.0),
            ("55_grid_refine", "level_1", 0.5),
        ]:
            src_dir = ctx.store.run_dir / "stage" / stage_name

            # Post-smooth mesh
            for suffix in [f"mesh_{level}.ply", f"mesh_{level}_ascii.ply"]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    # Also put in run root
                    shutil.copy2(src, run_root / suffix)

            # Pre-smooth mesh
            for suffix in [
                f"mesh_{level}_pre_smooth.ply",
                f"mesh_{level}_pre_smooth_ascii.ply",
            ]:
                src = src_dir / suffix
                if src.exists():
                    shutil.copy2(src, stage_dir / f"comparison_{suffix}")
                    shutil.copy2(src, run_root / suffix)

            if (src_dir / f"mesh_{level}.ply").exists():
                print(
                    f"[{self.key}]   copied {level} mesh ({voxel}A grid) + pre-smooth"
                )

```

ribctl/lib/npet2/__init__.py
```py

```

ribctl/lib/npet2/__main__.py
```py
# ribctl/lib/npet2/__main__.py
"""
npet2 CLI entry point.

Usage:
    python -m ribctl.lib.npet2 run 7K00 4UG0 --workers 4
    python -m ribctl.lib.npet2 run --from-file structures.txt --output-dir ./results
    python -m ribctl.lib.npet2 run 7K00 --cylinder-radius 40 --voxel-size 0.5
    python -m ribctl.lib.npet2 run 7K00 --mmcif /path/to/7K00.cif --profile /path/to/profile.json --landmarks /path/to/landmarks.json
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path
from typing import Optional


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="npet2",
        description="Ribosome exit tunnel geometry pipeline",
    )
    sub = p.add_subparsers(dest="command")

    # --- run ---
    run_p = sub.add_parser("run", help="Run the pipeline on one or more structures")

    # Input selection
    run_p.add_argument("rcsb_ids", nargs="*", help="RCSB IDs to process")
    run_p.add_argument("--from-file", type=str, default=None,
                       help="Read RCSB IDs from a text file (one per line)")

    # Provider mode
    run_p.add_argument("--mode", choices=["riboxyz", "standalone"], default="riboxyz",
                       help="riboxyz: use local riboxyz assets. standalone: use mmcif + profile files or API.")
    run_p.add_argument("--mmcif", type=str, default=None,
                       help="(standalone) Path to mmCIF file. For multiple structures, use a directory.")
    run_p.add_argument("--profile", type=str, default=None,
                       help="(standalone) Path to profile JSON file or directory of profiles.")
    run_p.add_argument("--landmarks", type=str, default=None,
                       help="(standalone) Path to landmarks JSON file or directory.")
    run_p.add_argument("--api-url", type=str, default=None,
                       help="(standalone) riboxyz API base URL for fetching profiles/landmarks.")

    # Output
    run_p.add_argument("--output-dir", type=str, default=None,
                       help="Root output directory (default: NPET2_RUNS_ROOT)")

    # Parallelism
    run_p.add_argument("--workers", "-j", type=int, default=1,
                       help="Number of parallel workers (default: 1)")

    # Config overrides (the commonly-tuned ones)
    cfg = run_p.add_argument_group("config overrides")
    cfg.add_argument("--cylinder-radius", type=float, default=None)
    cfg.add_argument("--cylinder-height", type=float, default=None)
    cfg.add_argument("--ptc-extension", type=float, default=None)
    cfg.add_argument("--voxel-size", type=float, default=None,
                     help="Level-0 grid voxel size in Angstroms")
    cfg.add_argument("--refine-voxel-size", type=float, default=None)
    cfg.add_argument("--no-mesh", action="store_true", help="Disable mesh generation")
    cfg.add_argument("--no-refine", action="store_true", help="Skip Stage55 grid refinement")
    cfg.add_argument("--dbscan-coarse-eps", type=float, default=None)
    cfg.add_argument("--dbscan-coarse-min-samples", type=int, default=None)
    cfg.add_argument("--dbscan-refine-eps", type=float, default=None)
    cfg.add_argument("--dbscan-refine-min-samples", type=int, default=None)
    cfg.add_argument("--config-json", type=str, default=None,
                     help="Path to a full RunConfig JSON (overrides individual flags)")

    # --- show-config ---
    sub.add_parser("show-config", help="Print the default RunConfig as JSON")

    return p


def _build_config(args) -> "RunConfig":
    from ribctl.lib.npet2.core.config import RunConfig, GridLevelConfig

    if args.config_json:
        import json
        data = json.loads(Path(args.config_json).read_text())
        # Reconstruct GridLevelConfig objects
        if "grid_levels" in data:
            data["grid_levels"] = [GridLevelConfig(**gl) for gl in data["grid_levels"]]
        return RunConfig(**data)

    kwargs = {}
    if args.cylinder_radius is not None:
        kwargs["cylinder_radius_A"] = args.cylinder_radius
    if args.cylinder_height is not None:
        kwargs["cylinder_height_A"] = args.cylinder_height
    if args.ptc_extension is not None:
        kwargs["cylinder_ptc_extension_A"] = args.ptc_extension
    if args.refine_voxel_size is not None:
        kwargs["refine_voxel_size_A"] = args.refine_voxel_size
    if args.no_mesh:
        kwargs["mesh_level0_enable"] = False
        kwargs["mesh_level1_enable"] = False
    if args.dbscan_coarse_eps is not None:
        kwargs["dbscan_level0_coarse_eps_A"] = args.dbscan_coarse_eps
    if args.dbscan_coarse_min_samples is not None:
        kwargs["dbscan_level0_coarse_min_samples"] = args.dbscan_coarse_min_samples
    if args.dbscan_refine_eps is not None:
        kwargs["dbscan_level0_refine_eps_A"] = args.dbscan_refine_eps
    if args.dbscan_refine_min_samples is not None:
        kwargs["dbscan_level0_refine_min_samples"] = args.dbscan_refine_min_samples

    if args.voxel_size is not None:
        kwargs["grid_levels"] = [
            GridLevelConfig(name="level_0", voxel_size_A=args.voxel_size,
                            occupancy_backend="legacy_kdtree"),
        ]

    return RunConfig(**kwargs)


def _collect_rcsb_ids(args) -> list[str]:
    ids = list(args.rcsb_ids) if args.rcsb_ids else []
    if args.from_file:
        p = Path(args.from_file)
        if not p.exists():
            print(f"Error: --from-file {p} does not exist", file=sys.stderr)
            sys.exit(1)
        for line in p.read_text().splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                ids.append(line)
    if not ids:
        print("Error: no RCSB IDs specified", file=sys.stderr)
        sys.exit(1)
    return [x.upper() for x in ids]


def _make_providers(args, rcsb_id: str):
    """Return (structure_provider, landmark_provider) for a given structure."""
    if args.mode == "riboxyz":
        from ribctl.lib.npet2.adapters.riboxyz_providers import (
            RiboxyzStructureProvider,
            RiboxyzLandmarkProvider,
        )
        return RiboxyzStructureProvider(), RiboxyzLandmarkProvider()

    # standalone mode
    from ribctl.lib.npet2.adapters.standalone_providers import (
        FileStructureProvider,
        FileLandmarkProvider,
    )

    api_base = args.api_url

    # Resolve mmcif path
    mmcif_path = None
    if args.mmcif:
        p = Path(args.mmcif)
        if p.is_dir():
            # Look for {RCSB_ID}.cif or {rcsb_id}.cif
            for candidate in [f"{rcsb_id}.cif", f"{rcsb_id.lower()}.cif"]:
                if (p / candidate).exists():
                    mmcif_path = p / candidate
                    break
            if mmcif_path is None:
                raise FileNotFoundError(
                    f"No mmCIF file found for {rcsb_id} in {p}. "
                    f"Expected {rcsb_id}.cif"
                )
        else:
            mmcif_path = p

    if mmcif_path is None:
        raise ValueError(f"--mmcif is required in standalone mode (for {rcsb_id})")

    # Resolve profile path
    profile_path = None
    if args.profile:
        p = Path(args.profile)
        if p.is_dir():
            for candidate in [f"{rcsb_id}_profile.json", f"{rcsb_id}.json"]:
                if (p / candidate).exists():
                    profile_path = p / candidate
                    break
        else:
            profile_path = p

    # Resolve landmarks path
    landmarks_path = None
    if args.landmarks:
        p = Path(args.landmarks)
        if p.is_dir():
            for candidate in [f"{rcsb_id}_landmarks.json", f"{rcsb_id}.json"]:
                if (p / candidate).exists():
                    landmarks_path = p / candidate
                    break
        else:
            landmarks_path = p

    return (
        FileStructureProvider(mmcif_path, profile_path=profile_path, api_base=api_base),
        FileLandmarkProvider(landmarks_path=landmarks_path, api_base=api_base),
    )


def _run_single(
    rcsb_id: str,
    args,
    config: "RunConfig",
    output_root: Optional[Path],
) -> dict:
    """Run pipeline for a single structure. Returns a result dict."""
    from ribctl.lib.npet2.run import run_npet2

    try:
        sp, lp = _make_providers(args, rcsb_id)

        # Allow output dir override
        if output_root:
            import ribctl.lib.npet2.core.settings as settings
            settings.NPET2_RUNS_ROOT = output_root

        ctx = run_npet2(rcsb_id, config, structure_provider=sp, landmark_provider=lp)
        return {
            "rcsb_id": rcsb_id,
            "status": "success",
            "run_dir": str(ctx.store.run_dir),
        }
    except Exception as e:
        return {
            "rcsb_id": rcsb_id,
            "status": "failed",
            "error": str(e),
            "traceback": traceback.format_exc(),
        }


def _run_worker(packed_args: tuple) -> dict:
    """Wrapper for ProcessPoolExecutor."""
    rcsb_id, args_ns, config_dict, output_root_str = packed_args
    from ribctl.lib.npet2.core.config import RunConfig, GridLevelConfig

    # Reconstruct config from dict
    if "grid_levels" in config_dict:
        config_dict["grid_levels"] = [GridLevelConfig(**gl) for gl in config_dict["grid_levels"]]
    config = RunConfig(**config_dict)

    output_root = Path(output_root_str) if output_root_str else None
    return _run_single(rcsb_id, args_ns, config, output_root)


def main():
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    if args.command == "show-config":
        from ribctl.lib.npet2.core.config import RunConfig
        cfg = RunConfig()
        print(json.dumps(asdict(cfg), indent=2))
        return

    if args.command == "run":
        rcsb_ids = _collect_rcsb_ids(args)
        config = _build_config(args)
        output_root = Path(args.output_dir) if args.output_dir else None

        n_workers = min(args.workers, len(rcsb_ids))

        if args.no_refine:
            # We need to modify the pipeline stages -- simplest way is a flag
            # that run.py checks. For now, store it in config as a workaround.
            pass

        print(f"npet2: processing {len(rcsb_ids)} structure(s), workers={n_workers}")

        if n_workers <= 1:
            results = []
            for rid in rcsb_ids:
                r = _run_single(rid, args, config, output_root)
                results.append(r)
                status = r["status"]
                print(f"  {rid}: {status}" + (f" -> {r.get('run_dir', '')}" if status == "success" else f" ({r.get('error', '')})"))
        else:
            # Serialize config for multiprocessing
            config_dict = asdict(config)
            output_root_str = str(output_root) if output_root else None

            packed = [
                (rid, args, config_dict, output_root_str)
                for rid in rcsb_ids
            ]

            results = []
            with ProcessPoolExecutor(max_workers=n_workers) as pool:
                futures = {pool.submit(_run_worker, p): p[0] for p in packed}
                for fut in as_completed(futures):
                    rid = futures[fut]
                    try:
                        r = fut.result()
                    except Exception as e:
                        r = {"rcsb_id": rid, "status": "failed", "error": str(e)}
                    results.append(r)
                    status = r["status"]
                    print(f"  {rid}: {status}" + (
                        f" -> {r.get('run_dir', '')}" if status == "success"
                        else f" ({r.get('error', '')})"
                    ))

        # Summary
        ok = sum(1 for r in results if r["status"] == "success")
        fail = len(results) - ok
        print(f"\nnpet2: {ok} succeeded, {fail} failed out of {len(results)}")

        if fail > 0:
            sys.exit(1)


if __name__ == "__main__":
    main()
```

ribctl/lib/npet2/run.py
```py
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
from typing import Optional

from ribctl.lib.npet2.adapters.riboxyz_providers import (
    RiboxyzLandmarkProvider,
    RiboxyzStructureProvider,
)
from ribctl.lib.npet2.core.config import RunConfig
from ribctl.lib.npet2.core.manifest import RunManifest
from ribctl.lib.npet2.core.run_id import compute_run_id
from ribctl.lib.npet2.core.settings import NPET2_ROOT, NPET2_RUNS_ROOT
from ribctl.lib.npet2.core.store import LocalRunStore
from ribctl.lib.npet2.core.types import StageContext
from ribctl.lib.npet2.core.pipeline import Pipeline


from ribctl.lib.npet2.stages.bootstrap import Stage00Inputs, Stage10Landmarks
from ribctl.lib.npet2.stages.grid_refine import Stage55GridRefine
from ribctl.lib.npet2.stages.legacy_minimal import (
    Stage20ExteriorShell,
    Stage30RegionAtoms,
    Stage40EmptySpace,
    Stage50Clustering,
    # Stage60SurfaceNormals,
    Stage70MeshValidate,
)


def _pipeline_version() -> str:
    return "npet2-dev"


def run_npet2(
    rcsb_id: str,
    config: Optional[RunConfig] = None,
    *,
    structure_provider=None,
    landmark_provider=None,
) -> StageContext:
    rcsb_id = rcsb_id.upper()
    config = config or RunConfig()


    structure_provider = structure_provider or RiboxyzStructureProvider()
    landmark_provider = landmark_provider or RiboxyzLandmarkProvider()

    config_resolved = asdict(config)
    inputs_fp = {
        "structure": structure_provider.fingerprint(rcsb_id),
        "landmarks": landmark_provider.fingerprint(rcsb_id),
    }

    struct_runs_dir = NPET2_RUNS_ROOT / rcsb_id
    struct_runs_dir.mkdir(parents=True, exist_ok=True)

    run_id = compute_run_id(
        rcsb_id=rcsb_id,
        pipeline_version=_pipeline_version(),
        inputs_fp=inputs_fp,
        config_resolved=config_resolved,
        runs_dir=struct_runs_dir,
    )

    run_dir = struct_runs_dir / run_id
    run_dir.mkdir(parents=True, exist_ok=True)

    manifest = RunManifest(
        rcsb_id=rcsb_id,
        run_id=run_id,
        pipeline_version=_pipeline_version(),
        inputs={"fingerprints": inputs_fp},
        config_resolved=config_resolved,
    )
    store = LocalRunStore(run_dir=run_dir, manifest=manifest)

    ctx = StageContext(
        run_id=run_id,
        rcsb_id=rcsb_id,
        config=config,
        store=store,
        inputs={
            "structure_provider": structure_provider,
            "landmark_provider": landmark_provider,
        },
    )

    # in run_npet2()
    from ribctl.lib.npet2.core.cache import LocalStageCache
    ctx.inputs["stage_cache"] = LocalStageCache(NPET2_ROOT / "cache")
    ctx.inputs["inputs_fp"] = inputs_fp

    pipeline = Pipeline(
        [
            Stage00Inputs(),
            Stage10Landmarks(),
            Stage20ExteriorShell(),
            Stage30RegionAtoms(),
            Stage40EmptySpace(),
            Stage50Clustering(),
            Stage55GridRefine(),   
            # Stage60SurfaceNormals(),
            Stage70MeshValidate(),
        ]
    )

    pipeline.run(ctx)
    return ctx

```


Actually it would be awesome if we could also have the `settings.py` and the `config.py` merged so all of the configuration is in a single place and is eventually configurable via some simple .env file in docker. Tell me if you see any other such simplifications that we can make.