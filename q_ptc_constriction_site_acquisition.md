Hey i have a pipeline that is supposed to detect the coordinates of the constriction site and the PTC in the strucutures of the ribosome, but it fails in a bunch of cases with these errors:
```
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  p3 ribctl/ribd.py  etl  get-all                                                                                   [npet_refactor]
2026-02-19T11:41:45.398003+0100 | DEBUG    | ribctl.lib.libtax:ensure_taxid_db_exists:32 - NCBI taxonomy database found at /Users/rtviii/dev/riboxyz/ncbi_taxonomy.sqlite
Usage: ribd.py etl get-all [OPTIONS]
Try 'ribd.py etl get-all --help' for help.

Error: Missing option '-t' / '--asset-type'. Choose from:
        MMCIF,
        STRUCTURE_PROFILE,
        CONSTRICTION_SITE,
        PTC,
        THUMBNAIL,
        NPET_MESH_ASCII,
        NPET_MESH,
        ALPHA_SHAPE,
        SRL,
        L7L12STALK,
        PSTALK
(venv) ᢹ saeta.rtviii[ dev/riboxyz ]  p3 ribctl/ribd.py  etl  get-all  -t CONSTRICTION_SITE -t PTC                                                      [npet_refactor]
2026-02-19T11:41:57.396875+0100 | DEBUG    | ribctl.lib.libtax:ensure_taxid_db_exists:32 - NCBI taxonomy database found at /Users/rtviii/dev/riboxyz/ncbi_taxonomy.sqlite
Found 2561 structures to process
Processing 257 chunks using 4 processes
2026-02-19T11:42:01.708429+0100 | INFO     | ribctl.asset_manager.asset_registry:wrapped:105 - Asset exists at /Users/rtviii/dev/RIBETL_DATA/5XY3/5XY3_CONSTRICTION_SITE.json, skipping
2026-02-19T11:42:01.708698+0100 | INFO     | ribctl.asset_manager.asset_registry:wrapped:105 - Asset exists at /Users/rtviii/dev/RIBETL_DATA/5XY3/5XY3_PTC.json, skipping
2026-02-19T11:42:01.710626+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 8JSH: Could not find uL4 or uL22 in 8JSH
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
    raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))
ValueError: Could not find uL4 or uL22 in 8JSH
2026-02-19T11:42:01.710626+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 8JSH: Could not find uL4 or uL22 in 8JSH
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
    raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))
ValueError: Could not find uL4 or uL22 in 8JSH
2026-02-19T11:42:01.712316+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 4DV0: Could not find uL4 or uL22 in 4DV0
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
    raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))
ValueError: Could not find uL4 or uL22 in 4DV0
2026-02-19T11:42:01.712316+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 4DV0: Could not find uL4 or uL22 in 4DV0
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
    raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))
ValueError: Could not find uL4 or uL22 in 4DV0
2026-02-19T11:42:01.755844+0100 | INFO     | ribctl.asset_manager.asset_registry:wrapped:105 - Asset exists at /Users/rtviii/dev/RIBETL_DATA/5Z3G/5Z3G_PTC.json, skipping
2026-02-19T11:42:01.756188+0100 | INFO     | ribctl.asset_manager.asset_registry:wrapped:105 - Asset exists at /Users/rtviii/dev/RIBETL_DATA/5Z3G/5Z3G_CONSTRICTION_SITE.json, skipping
2026-02-19T11:42:04.395240+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_ptc for 6WDR: No LSU rRNA found in structure
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 158, in generate_ptc
    return PTC_location(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/ptc_via_trna.py", line 164, in PTC_location
    LSU_RNA_tgt_aaid = RO.get_LSU_rRNA().auth_asym_id
  File "/Users/rtviii/dev/riboxyz/ribctl/ribosome_ops.py", line 153, in get_LSU_rRNA
    raise Exception("No LSU rRNA found in structure")
Exception: No LSU rRNA found in structure
2026-02-19T11:42:04.395240+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_ptc for 6WDR: No LSU rRNA found in structure
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 158, in generate_ptc
    return PTC_location(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/ptc_via_trna.py", line 164, in PTC_location
    LSU_RNA_tgt_aaid = RO.get_LSU_rRNA().auth_asym_id
  File "/Users/rtviii/dev/riboxyz/ribctl/ribosome_ops.py", line 153, in get_LSU_rRNA
    raise Exception("No LSU rRNA found in structure")
Exception: No LSU rRNA found in structure
2026-02-19T11:42:04.399714+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 6WDR: Could not find uL4 or uL22 in 6WDR
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
    raise ValueError("Could not find uL4 or uL22 in {}".format(rcsb_id))
ValueError: Could not find uL4 or uL22 in 6WDR
2026-02-19T11:42:04.399714+0100 | ERROR    | ribctl.asset_manager.asset_registry:wrapped:114 - Failed generate_constriction for 6WDR: Could not find uL4 or uL22 in 6WDR
Traceback (most recent call last):
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 109, in wrapped
    result = await func(rcsb_id)
  File "/Users/rtviii/dev/riboxyz/ribctl/asset_manager/asset_registry.py", line 163, in generate_constriction
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())
  File "/Users/rtviii/dev/riboxyz/ribctl/lib/landmarks/constriction_site.py", line 17, in get_constriction
```
 
Can we debug?


Here is my repo:
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
├── npet2.env
├── pipeline_manager.py
├── PLAN_refactor_npet_pipeline_6.md
├── PLAN_refactor_npet_pipeline_7_cleanup.md
├── PLAN_refactor_npet_pipeline_8.md
├── q_mass_runs_and_packaging.md
├── q_ptc_constriction_site_acquisition.md
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
│   │   │   │   ├── geometry.py
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
e  [error opening dir]

35 directories, 182 files
```


Feel free to ask me for more files, but here are the main ones:
ribctl/asset_manager/asset_registry.py
```py
from typing import TypeVar, Callable, Awaitable
import functools
from loguru import logger
from pydantic import BaseModel
from pathlib import Path
from loguru import logger
from typing import Dict, Callable, Awaitable
from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.lib.npet.alphalib import alpha_contour_via_poisson_recon
# from ribctl.lib.npet.npet_driver import create_npet_mesh
from ribctl.lib.npet.npet_pipeline import create_npet_mesh
from ribctl.lib.utils import download_unpack_place
from ribctl.asset_manager.asset_types import AssetType
from ribctl import RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.etl.etl_collector import ETLCollector
from ribctl.lib.landmarks.constriction_site import get_constriction
from ribctl.lib.landmarks.ptc_via_trna import PTC_location
from ribctl.lib.schema.types_ribosome import (
    ConstrictionSite,
    PTCInfo,
    RibosomeStructure,
)

from .asset_types import AssetType
ModelT = TypeVar("ModelT", bound=BaseModel)

class RawAssetHandler:
    """Handler for raw file assets with extensible asset type matching"""

    def __init__(self):
        self._handlers: Dict[AssetType, Callable[[str, bool], Awaitable[None]]] = {}
        self._register_default_handlers()

    def _register_default_handlers(self) -> None:
        """Register built-in handlers for known raw asset types"""
        self.register_handler(AssetType.MMCIF, self._fetch_mmcif)
        self.register_handler(AssetType.NPET_MESH, npet_mesh_handler)
        self.register_handler(AssetType.ALPHA_SHAPE, alphashape_handler)

    def register_handler(
        self, asset_type: AssetType, handler: Callable[[str, bool], Awaitable[None]]
    ) -> None:
        """Register a new handler for an asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(
                f"Cannot register handler for non-raw asset type: {asset_type}"
            )
        self._handlers[asset_type] = handler

    async def handle_asset(
        self, rcsb_id: str, asset_type: AssetType, force: bool = False
    ) -> None:
        """Generic handler for any registered raw asset type"""
        if not asset_type.is_raw_asset:
            raise ValueError(f"Asset type {asset_type} is not a raw asset")

        handler = self._handlers.get(asset_type)
        if not handler:
            raise ValueError(f"No handler registered for raw asset type: {asset_type}")

        await handler(rcsb_id, force)

    async def _fetch_mmcif(self, rcsb_id: str, force: bool = False) -> None:
        """Download and save mmCIF file"""
        output_path = AssetType.MMCIF.get_path(rcsb_id)

        if output_path.exists() and not force:
            logger.info(f"MMCIF exists for {rcsb_id}, skipping")
            return

        output_path.parent.mkdir(parents=True, exist_ok=True)
        await download_unpack_place(rcsb_id)
        logger.success(f"Downloaded MMCIF for {rcsb_id}")

    # Example of how to add another handler:
    # async def _fetch_npet_mesh(self, rcsb_id: str, force: bool = False) -> None:
    #     """Download and save NPET mesh file"""
    #     output_path = self.base_dir / rcsb_id.upper() / "TUNNELS" / f"{rcsb_id}_NPET_MESH.ply"
    #
    #     if output_path.exists() and not force:
    #         logger.info(f"NPET mesh exists for {rcsb_id}, skipping")
    #         return
    #
    #     output_path.parent.mkdir(parents=True, exist_ok=True)
    #     # Add actual download/generation logic here
    #     logger.success(f"Generated NPET mesh for {rcsb_id}")


class AssetRegistry:
    def __init__(self, manager: RibosomeAssetManager):
        self.manager = manager
        self.raw_handler = RawAssetHandler()

    def register(self, asset_type: AssetType):
        def decorator(
            func: Callable[[str], Awaitable[ModelT]]
        ) -> Callable[[str, bool], Awaitable[None]]:
            @functools.wraps(func)
            async def wrapped(rcsb_id: str, overwrite: bool = False) -> None:
                output_path = asset_type.get_path(rcsb_id)
                try:
                    if output_path.exists() and not overwrite:
                        logger.info(f"Asset exists at {output_path}, skipping")
                        return

                    output_path.parent.mkdir(parents=True, exist_ok=True)
                    result = await func(rcsb_id)
                    output_path.write_text(result.model_dump_json())
                    logger.success(f"Generated {asset_type.name} for {rcsb_id}")

                except Exception as e:
                    logger.exception(f"Failed {func.__name__} for {rcsb_id}: {str(e)}")
                    raise

            self.manager.register_generator(asset_type, wrapped)
            return wrapped

        return decorator

    async def generate_asset(
        self, rcsb_id: str, asset_type: AssetType, force: bool = False
    ) -> None:
        if asset_type.is_raw_asset:
            await self.raw_handler.handle_asset(rcsb_id, asset_type, force)
        else:
            asset_def = self.manager.assets[asset_type]
            if not asset_def.generator:
                raise ValueError(f"No generator registered for {asset_type}")
            await asset_def.generator(rcsb_id, force)

    async def generate_multiple(
        self, rcsb_id: str, asset_types: list[AssetType], force: bool = False
    ) -> None:
        """Generate multiple assets for a structure"""
        for asset_type in asset_types:
            await self.generate_asset(rcsb_id, asset_type, force)

async def npet_mesh_handler(rcsb_id: str, force: bool) -> None:
    create_npet_mesh(rcsb_id, Path('/Users/rtviii/dev/riboxyz/ribctl/lib/npet/pipeline/logs'))

async def alphashape_handler(rcsb_id: str, force: bool) -> None:
    alpha_contour_via_poisson_recon(rcsb_id)

main_registry = AssetRegistry(RibosomeAssetManager(RIBETL_DATA))

@main_registry.register(AssetType.STRUCTURE_PROFILE)
async def generate_profile(rcsb_id: str) -> RibosomeStructure:
    profile = await ETLCollector(rcsb_id).generate_profile(
        overwrite=False, reclassify=True
    )
    return profile


@main_registry.register(AssetType.PTC)
async def generate_ptc(rcsb_id: str) -> PTCInfo:
    return PTC_location(rcsb_id)


@main_registry.register(AssetType.CONSTRICTION_SITE)
async def generate_constriction(rcsb_id: str) -> ConstrictionSite:
    return ConstrictionSite(location=get_constriction(rcsb_id).tolist())

```

ribctl/asset_manager/asset_manager.py
```py
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

```

ribctl/asset_manager/asset_types.py
```py
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
        return Path(
            self.value.path_template.format(asset_dir=asset_dir, pdb_id=pdb_id.upper())
        )

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

    NPET_MESH_ASCII = AssetInfo(
        "npet_mesh_ascii",
        path_template="{asset_dir}/{pdb_id}_NPET_MESH_ascii.ply",
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

    NPET_MESH = AssetInfo(
        "npet_mesh",
        path_template="{asset_dir}/{pdb_id}_NPET_MESH.ply",
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
        is_raw=True,
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

```

ribctl/asset_manager/assets_structure.py
```py
from enum import   auto
import json
import os
from pprint import pprint
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Structure import Structure
from Bio.PDB.MMCIFParser import FastMMCIFParser
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass,  RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA

class StructureAssetPaths:
    rcsb_id:str
    def __init__(self, rcsb_id) -> None:
        self.rcsb_id = rcsb_id
        pass

    def nonpoly_entity(self, chemId:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_{chemId.upper()}_STRUCTURE.cif"

    def binding_site(self, chemId:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}.json"

    def binding_site_prediction(self, chemId:str, source_struct:str):
        return f"{self.dir}/{self.rcsb_id.upper()}_LIG_{chemId.upper()}_PREDICTION_VIA_{source_struct.upper()}.json"
    
    @property
    def cif(self):
        return os.path.join(RIBETL_DATA, self.rcsb_id, f"{self.rcsb_id}.cif")

    @property
    def dir(self):
        return os.path.join(RIBETL_DATA, self.rcsb_id)

    @property
    def profile(self):
        return os.path.join(self.dir, f"{self.rcsb_id}.json")

    @property
    def chains_dir(self):
        return f"{self.dir}/CHAINS"

    @property
    def classification_report(self):
        return os.path.join(self.dir, f"classification_report_{self.rcsb_id}.json")

    @property
    def thumbnail(self):
        return f"{self.dir}/{self.rcsb_id}.png"

class StructureAssets:

    rcsb_id: str
    paths  : StructureAssetPaths

    def __init__(self, rcsb_id: str) -> None:
        self.rcsb_id = rcsb_id
        self.paths   = StructureAssetPaths(rcsb_id)
    
    def biopython_structure(self)-> Structure:
        cifpath = StructureAssetPaths(self.rcsb_id).cif
        return FastMMCIFParser(QUIET=True).get_structure(self.rcsb_id, cifpath)

    def profile(self) -> RibosomeStructure:
        with open(self.paths.profile, "r") as f:
            return RibosomeStructure.model_validate(json.load(f))

    def biopython_get_chain(self, auth_asym_id: str) -> Chain:
        return self.biopython_structure().child_dict[0].child_dict[auth_asym_id]

    def _verify_dir_exists(self):
        if not os.path.exists(self.paths.dir):
            os.umask(0)
            os.makedirs(self.paths.dir, 0o755)
```

ribctl/ribd.py
```py
import os
from pprint import pprint
import sys
from loguru import logger
from ribctl.lib.libtax import ensure_taxid_db_exists

sys.dont_write_bytecode = True
sys.path.append("/home/rtviii/dev/riboxyz")
from ribctl.lib.libbsite import extract_ligand_to_mmcif
from logger_config import configure_logging

configure_logging()
import concurrent.futures as cf
from neo4j_ribosome.db_lib_reader import Neo4jReader
from functools import partial
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_registry import AssetRegistry
from ribctl.lib.schema.types_ribosome import PTCInfo, RibosomeStructure
from ribctl.asset_manager.assets_structure import StructureAssets
from ribctl.global_ops import GlobalOps
from ribctl.ribosome_ops import RibosomeOps
from neo4j_ribosome import NEO4J_CURRENTDB, NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER
from neo4j_ribosome.db_lib_builder import Neo4jAdapter
from ribctl import RIBETL_DATA
from typing import List, Optional, Tuple, Dict
import click
from pathlib import Path
import asyncio
from concurrent.futures import ProcessPoolExecutor
import math
from ribctl.asset_manager.parallel_acquisition import (
    process_chunk,
    AcquisitionResult,
    process_chunk_with_tracking,
)
import multiprocessing
from ribctl.asset_manager.asset_registry import main_registry
from concurrent.futures import (
    ProcessPoolExecutor,
)
from tqdm import tqdm
from ribctl.asset_manager.asset_types import AssetType


def get_input_pdb_ids() -> List[str]:
    """Get PDB IDs from either stdin (if piped) or return None to handle as argument"""
    if not sys.stdin.isatty():
        return [line.strip().upper() for line in sys.stdin if line.strip()]
    return []


class PDBIDsParam(click.ParamType):
    """Custom parameter type for PDB IDs validation"""

    name = "pdb_ids"

    def convert(self, value, param, ctx):

        if type(value) == list:  # Changed from isinstance(value, list)
            return value
        if not value:
            return []
        pdb_id = value.strip().upper()
        if not (len(pdb_id) == 4 and pdb_id.isalnum()):
            self.fail(f"Invalid PDB ID format: {value}", param, ctx)
        return [pdb_id]


def get_input_structures() -> List[str]:
    """Read PDB IDs from stdin if available, otherwise return empty list"""
    if not sys.stdin.isatty():
        return [line.strip().upper() for line in sys.stdin.readlines() if line.strip()]
    return []


@click.group()
@click.pass_context
def cli(ctx):
    """ribctl - Command line interface for ribosome data pipeline"""
    ctx.ensure_object(dict)
    ctx.obj["piped_pdb_ids"] = get_input_pdb_ids()

@cli.group()
@click.pass_context
def etl(ctx):
    """ETL operations for ribosome data"""
    # Run the check here so every ETL command is safe
    ensure_taxid_db_exists()
    pass

@etl.command()
@click.pass_context
@click.option(
    "-t",
    "--asset-type",
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    multiple=True,
    required=True,
    help="Asset type(s) to acquire. Can be specified multiple times.",
)
@click.argument("pdb_id", required=False)
@click.option("--force", is_flag=True, help="Force regeneration of existing assets")
@click.option(
    "--max-structures",
    default=4,
    help="Maximum number of concurrent structures to process",
)
@click.option(
    "--max-assets", default=3, help="Maximum number of concurrent assets per structure"
)
def get(
    ctx,
    asset_type: List[str],
    pdb_id: str,
    force: bool,
    max_structures: int,
    max_assets: int,
):
    """
    Acquire assets for specified structure(s). Accepts either a single PDB ID as argument
    or a list of IDs via stdin pipe.
    """
    # Get structures from either argument or pipe
    structures = ctx.obj.get("piped_pdb_ids", [])
    
    # If nothing in context, try reading from stdin directly
    if not structures:
        structures = get_input_structures()
        
    if pdb_id:
        structures.append(pdb_id.upper())
    

    # Convert asset type names to enum
    asset_types = [AssetType[t] for t in asset_type]

    async def process_structure(rcsb_id: str):
        try:
            await main_registry.generate_multiple(rcsb_id, asset_types, force)
            click.echo(f"Successfully processed assets for {rcsb_id}")
        except Exception as e:
            logger.exception(f"Failed to process {rcsb_id}")
            click.echo(f"Failed to process {rcsb_id}: {str(e)}", err=True)

    async def process_all():
        sem = asyncio.Semaphore(max_structures)

        async def wrapped_process(rcsb_id: str):
            async with sem:
                await process_structure(rcsb_id)

        tasks = [wrapped_process(rcsb_id) for rcsb_id in structures]
        await asyncio.gather(*tasks)

    # Run the async processing
    asyncio.run(process_all())


@etl.command()
@click.option(
    "-t",
    "--asset-type",
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    multiple=True,
    required=True,
    help="Asset type(s) to acquire. Can be specified multiple times.",
)
@click.option("--force", is_flag=True, help="Force regeneration of existing assets")
@click.option(
    "--max-structures",
    default=4,
    help="Maximum number of concurrent structures per process",
)
@click.option(
    "--max-assets", default=3, help="Maximum number of concurrent assets per structure"
)
@click.option("--processes", default=4, help="Number of parallel processes to use")
@click.option(
    "--chunk-size", default=10, help="Number of structures to process per chunk"
)
def get_all(
    asset_type: List[str],
    force: bool,
    max_structures: int,
    max_assets: int,
    processes: int,
    chunk_size: int,
):
    """
    Acquire assets for ALL structures in the database using parallel processing.
    Structures are processed in parallel chunks using multiple CPU cores.
    """
    # Convert asset type names to enum
    asset_types = [AssetType[t] for t in asset_type]

    # Get all available structures
    structures = GlobalOps.list_profiles()

    if not structures:
        click.echo("No structures found in the database!", err=True)
        return

    total_structures = len(structures)
    click.echo(f"Found {total_structures} structures to process")

    # Split structures into chunks
    chunks = [
        structures[i : i + chunk_size] for i in range(0, len(structures), chunk_size)
    ]

    # Create a partial function with all the fixed arguments
    process_func = partial(
        process_chunk_with_tracking,
        base_dir=str(Path(RIBETL_DATA)),
        asset_type_names=[t.name for t in asset_types],
        force=force,
        max_structures=max_structures,
        max_assets=max_assets,
    )

    # Process chunks in parallel
    click.echo(f"Processing {len(chunks)} chunks using {processes} processes")
    processed = 0
    failed = 0

    with ProcessPoolExecutor(max_workers=processes) as executor:
        for result in executor.map(process_func, chunks):
            # Count successes and failures from results
            for rcsb_id, acquisitions in result.items():
                if all(acq.success for acq in acquisitions):
                    processed += 1
                else:
                    failed += 1

    # Print final summary
    click.echo(f"\nProcessing complete:")
    click.echo(f"Total structures: {total_structures}")
    click.echo(f"Successfully processed: {processed}")
    click.echo(f"Failed: {failed}")

    # Return non-zero exit code if any failures
    if failed > 0:
        raise click.ClickException(f"Failed to process {failed} structures")


@etl.command()
@click.argument("pdb_ids", type=PDBIDsParam(), nargs=-1)
@click.pass_context
def verify(ctx, pdb_ids):
    """Verify assets exist for given PDB IDs"""
    all_pdb_ids = list(set(ctx.obj["piped_pdb_ids"] + sum(pdb_ids, [])))
    if not all_pdb_ids:
        click.echo("No PDB IDs provided", err=True)
        return




@etl.command()
@click.pass_context
def verify_all(ctx):
    """Verify all assets in the system by validating each structure against the RibosomeStructure model"""
    from concurrent.futures import ThreadPoolExecutor
    import asyncio
    from typing import Dict, List, Tuple

    structures = GlobalOps.list_profiles()
    if not structures:
        click.echo("No structures found in the system", err=True)
        return

    click.echo(f"Found {len(structures)} structures to verify")

    # Track validation results
    results: Dict[str, List[Tuple[str, str]]] = {
        "valid": [],
        "invalid": []
    }

    async def validate_structure(rcsb_id: str) -> None:
        """Validate a single structure and update results"""
        try:
            ops = RibosomeOps(rcsb_id)
            structure_data = ops.assets.profile()
            
            # Attempt validation
            RibosomeStructure.model_validate(structure_data)
            results["valid"].append(rcsb_id)
            
        except Exception as e:
            results["invalid"].append((rcsb_id, str(e)))

    async def process_all_structures():
        """Process all structures concurrently with a progress bar"""
        tasks = []
        with click.progressbar(
            structures,
            label="Validating structures",
            length=len(structures)
        ) as progress_structures:
            for rcsb_id in progress_structures:
                task = asyncio.create_task(validate_structure(rcsb_id))
                tasks.append(task)
            await asyncio.gather(*tasks)

    # Run validation
    try:
        asyncio.run(process_all_structures())
    except Exception as e:
        click.echo(f"\nError during validation: {str(e)}", err=True)
        return

    # Report results
    total = len(structures)
    valid_count = len(results["valid"])
    invalid_count = len(results["invalid"])

    click.echo("\nValidation Summary:")
    click.echo(f"Total structures processed: {total}")
    click.echo(f"Valid structures: {valid_count}")
    click.echo(f"Invalid structures: {invalid_count}")

    if invalid_count > 0:
        click.echo("\nInvalid structures:")
        for rcsb_id, error in results["invalid"]:
            click.echo(f"  - {rcsb_id}: {error}")

    # Return non-zero exit code if any structures failed validation
    if invalid_count > 0:
        ctx.exit(1)



@etl.command()
@click.pass_context
def ligand_mmcifs(ctx):
    from neo4j_ribosome.db_lib_reader import dbqueries

    ligs = dbqueries.ligands_per_structure()
    for lig in ligs:
        chemid = lig["chemicalId"]
        structs = lig["structures"]
        for struct in structs:
            structure = RibosomeOps(struct).assets.biopython_structure()
            try:
                if os.path.exists(StructureAssets(struct).paths.nonpoly_entity(chemid)):
                    print("Already exists: ", struct, chemid)
                    continue
                extract_ligand_to_mmcif(
                    structure,
                    chemid,
                    StructureAssets(struct).paths.nonpoly_entity(chemid),
                )
            except:
                pass
# Define this outside the command function to make it picklable
# For multiprocessing
import concurrent.futures
from tqdm import tqdm
import sys
import time

def validate_structure(rcsb_id):
    """Validate a single structure against the Pydantic model"""
    from pydantic import ValidationError
    
    try:
        ops = RibosomeOps(rcsb_id)
        # Get structure profile
        structure_data = ops.assets.profile()
        
        # Validate against Pydantic model
        RibosomeStructure.model_validate(structure_data)
        
        return {
            "rcsb_id": rcsb_id,
            "valid": True,
            "error": None
        }
    except ValidationError as e:
        return {
            "rcsb_id": rcsb_id,
            "valid": False,
            "error": str(e)
        }
    except Exception as e:
        return {
            "rcsb_id": rcsb_id,
            "valid": False,
            "error": f"Failed to load structure: {str(e)}"
        }


@etl.command()
@click.option(
    "--show-errors",
    is_flag=True,
    help="Show detailed Pydantic validation errors for failed structures",
)
@click.option(
    "--verbose", "-v",
    is_flag=True,
    help="Show detailed progress information",
)
@click.option(
    "--workers",
    default=max(1, multiprocessing.cpu_count() - 1),
    help="Number of worker processes for parallel validation",
)
def check_schema(show_errors, verbose, workers):
    """
    Validate all structure profiles against the Pydantic model schema.
    
    Attempts to load each structure profile and validates it against the
    RibosomeStructure Pydantic model. Reports structures that fail validation.
    
    Examples:\n
    \b
    # Check all structures with minimal output
    ribd.py etl check_schema
    
    \b
    # Show detailed validation errors
    ribd.py etl check_schema --show-errors
    
    \b
    # Verbose output with 8 worker processes
    ribd.py etl check_schema -v --workers 8
    """
    from concurrent.futures import ProcessPoolExecutor
    import time
    import concurrent.futures
    
    # Get all available structures
    structures = GlobalOps.list_profiles()
    total = len(structures)
    
    if not structures:
        click.echo("No structures found.")
        return
    
    click.echo(f"Validating schema for {total} structures...")
    
    # Results tracking
    results = {
        "passed": [],
        "failed": [],
        "errors": {}
    }
    
    start_time = time.time()
    
    # Process structures in parallel
    with ProcessPoolExecutor(max_workers=workers) as executor:
        with tqdm(total=total, disable=not verbose) as progress_bar:
            futures = [executor.submit(validate_structure, rcsb_id) for rcsb_id in structures]
            
            for future in concurrent.futures.as_completed(futures):
                result = future.result()
                rcsb_id = result["rcsb_id"]
                
                if result["valid"]:
                    results["passed"].append(rcsb_id)
                    if verbose:
                        click.echo(f"✓ {rcsb_id}: Valid")
                else:
                    results["failed"].append(rcsb_id)
                    results["errors"][rcsb_id] = result["error"]
                    if verbose:
                        click.echo(f"✗ {rcsb_id}: Invalid")
                
                progress_bar.update(1)
    
    # Print summary
    elapsed_time = time.time() - start_time
    passed_count = len(results["passed"])
    failed_count = len(results["failed"])
    
    click.echo("\nValidation Summary:")
    click.echo(f"Total structures: {total}")
    click.echo(f"Passed: {passed_count} ({passed_count/total*100:.1f}%)")
    click.echo(f"Failed: {failed_count} ({failed_count/total*100:.1f}%)")
    click.echo(f"Time elapsed: {elapsed_time:.2f} seconds")
    
    # Print failed structures
    if failed_count > 0:
        click.echo("\nFailed structures:")
        for i, rcsb_id in enumerate(sorted(results["failed"]), 1):
            click.echo(f"{rcsb_id}")
            
            # Print detailed errors if requested
            if show_errors:
                error_msg = results["errors"][rcsb_id]
                formatted_error = "\n".join(f"    {line}" for line in error_msg.split("\n"))
                click.echo(f"   Error details:\n{formatted_error}\n")
    
    # Return non-zero exit code if any structures failed validation
    if failed_count > 0:
        sys.exit(1)

@etl.command(name="sync_all")
@click.option(
    "--asset-types",
    "-t",
    multiple=True,
    type=click.Choice([t.name for t in AssetType], case_sensitive=True),
    help="Asset types to acquire (required)",
    required=True,
)
@click.option(
    "--force", "-f", is_flag=True, help="Force regeneration of existing assets"
)
@click.option(
    "--workers",
    "-w",
    default=max(1, multiprocessing.cpu_count() - 1),
    help="Number of worker processes",
)
@click.option(
    "--chunk-size", "-c", default=4, help="Number of structures to process per worker"
)
@click.option(
    "--concurrent-structures",
    "-s",
    default=4,
    help="Maximum concurrent structures per worker",
)
@click.option(
    "--concurrent-assets",
    "-a",
    default=3,
    help="Maximum concurrent assets per structure",
)
@click.option("--delay", "-d", default=2.0, help="Delay in seconds between chunks")
@click.option(
    "--max-retries",
    "-r",
    default=3,
    help="Maximum number of retries per failed operation",
)
@click.pass_context
def sync_all(
    ctx,
    asset_types,
    force,
    workers,
    chunk_size,
    concurrent_structures,
    concurrent_assets,
    delay,
    max_retries,
):
    """Synchronize all structure assets with RCSB, downloading or updating as needed.

    Compares local assets with RCSB database and only processes structures that:
    - Don't exist locally
    - Have different modification dates
    - Are missing requested assets

    Examples:\n
    \b
    # Sync specific assets for all outdated structures
    ribd.py etl sync_all --asset-types MMCIF STRUCTURE_PROFILE

    \b
    # Force sync with custom settings
    ribd.py etl sync_all --asset-types MMCIF PTC --workers 4 --force
    """
    from time import sleep
    import random

    ensure_taxid_db_exists()
    def process_with_retry(operation, max_retries):
        """Execute operation with retries and exponential backoff"""
        for attempt in range(max_retries):
            try:
                return operation()
            except Exception as e:
                if attempt == max_retries - 1:
                    raise e
                sleep_time = (attempt + 1) * 2 + random.random()
                sleep(sleep_time)

    try:
        # Get structures that need updating
        click.echo("Comparing local assets with RCSB database...")

        # Here we'd use your vs_rcsb method to get structures needing updates
        structures_to_update = GlobalOps.missing_profiles()

        if not structures_to_update:
            click.echo("All structures are up to date!")
            return

        click.echo(f"Found {len(structures_to_update)} structures needing updates")
    except Exception as e:
        click.echo(f"Failed to compare with RCSB: {str(e)}", err=True)
        return

    # Convert asset type names to enums
    selected_asset_types = [AssetType[name] for name in asset_types]
    assets_str = ", ".join(ast.name for ast in selected_asset_types)

    # Inform user about the operation
    click.echo(
        f"Getting assets [{assets_str}] for {len(structures_to_update)} structures..."
    )
    click.echo(f"Using {workers} workers, {chunk_size} structures per chunk")
    click.echo(f"Delay between chunks: {delay}s")

    # Split structures into chunks
    chunks = [
        structures_to_update[i : i + chunk_size]
        for i in range(0, len(structures_to_update), chunk_size)
    ]

    errors = []
    with tqdm(total=len(structures_to_update)) as pbar:
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = []

            # Submit all chunks to the process pool
            for chunk in chunks:
                future = executor.submit(
                    process_chunk,
                    str(RIBETL_DATA),
                    chunk,
                    [ast.name for ast in selected_asset_types],
                    force,
                    concurrent_structures,
                    concurrent_assets,
                )
                futures.append(future)

            # Process results as they complete
            for idx, future in enumerate(futures):
                try:
                    results = process_with_retry(lambda: future.result(), max_retries)

                    # Process results for this chunk
                    for rcsb_id, asset_results in results.items():
                        pbar.update(1)
                        pbar.set_description(f"Processed {rcsb_id}")

                        # Report failures
                        for result in asset_results:
                            if not result.success:
                                error_msg = (
                                    f"Error with {result.asset_type_name} "
                                    f"for {rcsb_id}: {result.error}"
                                )
                                errors.append((rcsb_id, error_msg))
                                click.echo(f"\n{error_msg}", err=True)

                    # Apply delay after each chunk except the last
                    if idx < len(futures) - 1:
                        sleep(delay)

                except Exception as e:
                    error_msg = (
                        f"Failed to process chunk after {max_retries} retries: {str(e)}"
                    )
                    click.echo(f"\n{error_msg}", err=True)
                    chunk = chunks[idx]  # Get corresponding chunk for this future
                    errors.extend((str(pdb_id), error_msg) for pdb_id in chunk)

    # Final report
    total_processed = len(structures_to_update)
    failed = len(errors)
    succeeded = total_processed - failed

    click.echo("\nSync Summary:")
    click.echo(f"Total structures processed: {total_processed}")
    click.echo(f"Successfully processed: {succeeded}")
    click.echo(f"Failed: {failed}")

    if errors:
        click.echo("\nErrors occurred during processing:")
        for rcsb_id, error in errors:
            click.echo(f"  - {rcsb_id}: {error}")
    else:
        click.echo("\nAll structures processed successfully")


@cli.group()
@click.pass_context
def db(ctx):
    """Database operations for ribosome data"""
    pass


@db.command()
@click.argument("query_file", type=click.Path(exists=True))
@click.option(
    "--params",
    "-p",
    multiple=True,
    help="Parameters for the Cypher query in key=value format",
)
@click.option("--output", "-o", type=click.Path(), help="Output file for query results")
@click.pass_context
def cypher(ctx, query_file, params, output):
    """Execute a Cypher query from a file"""
    # Parse parameters
    query_params = {}
    for param in params:
        try:
            key, value = param.split("=")
            query_params[key.strip()] = value.strip()
        except ValueError:
            click.echo(
                f"Invalid parameter format: {param}. Use key=value format.", err=True
            )
            return

    # Read query from file
    with open(query_file, "r") as f:
        query = f.read()

    click.echo(f"Executing Cypher query from {query_file}")
    # TODO: Implement query execution logic
    # results = execute_cypher_query(query, query_params)

    # if output:
    #     # Save results to file
    #     with open(output, 'w') as f:
    #         json.dump(results, f)
    # else:
    #     # Print results to stdout
    #     click.echo(results)


@db.command()
@click.option("--force", "-f", is_flag=True, help="Force reinitialization of database")
@click.argument("instance_name", required=True)
@click.pass_context
def init(ctx, force, instance_name):
    """Initialize a new Neo4j database instance.

    Creates necessary constraints and initial data structures in the specified instance.

    Examples:\n
    \b
    # Initialize a new instance named 'ribosome'
    ribd.py db init ribosome

    \b
    # Force reinitialize an existing instance
    ribd.py db init --force ribosome
    """
    if not force:
        click.confirm(
            f'This will initialize database instance "{instance_name}". Continue?',
            abort=True,
        )

    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, instance_name, NEO4J_PASSWORD)
        adapter.initialize_new_instance()
        click.echo(f"Database instance '{instance_name}' initialized successfully")
    except Exception as e:
        click.echo(f"Failed to initialize database instance: {str(e)}", err=True)


@db.command()
@click.argument("pdb_id", type=str)
@click.option(
    "--force", is_flag=True, help="Force upload even if structure already exists"
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload(
    pdb_id: str,
    force: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload a structure to the database using its PDB ID.

    PDB_ID should be a valid RCSB PDB identifier (e.g., '1j5e').
    """
    try:
        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Check if structure exists (unless force flag is used)
        if not force and adapter.check_structure_exists(pdb_id):
            click.echo(
                f"Structure {pdb_id.upper()} already exists in database. Use --force to overwrite."
            )
            return

        # Upload the structure
        click.echo(f"Uploading structure {pdb_id.upper()}...")
        adapter.add_total_structure(pdb_id, disable_exists_check=force)
        click.echo(f"Successfully uploaded structure {pdb_id.upper()}")

    except Exception as e:
        click.echo(f"Error uploading structure {pdb_id.upper()}: {str(e)}", err=True)
        raise click.Abort()


@db.command()
@click.option("--workers", "-w", type=int, default=4, help="Number of parallel workers")
@click.option(
    "--force", is_flag=True, help="Force upload even if structures already exist"
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload_all(
    workers: int,
    force: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload all available structures to the database in parallel.

    Uses GlobalOps to get list of structures and uploads them using multiple workers.
    """

    async def upload_structure(
        adapter: Neo4jAdapter, pdb_id: str, force: bool
    ) -> tuple[str, bool, Optional[str]]:
        """Upload a single structure and return result"""
        try:
            if not force and adapter.check_structure_exists(pdb_id):
                return pdb_id, False, "Already exists"

            adapter.add_total_structure(pdb_id, disable_exists_check=force)
            return pdb_id, True, None

        except Exception as e:
            return pdb_id, False, str(e)

    async def upload_worker(
        queue: asyncio.Queue, adapter: Neo4jAdapter, force: bool
    ) -> None:
        """Worker to process structures from the queue"""
        while True:
            try:
                pdb_id = await queue.get()
                result = await upload_structure(adapter, pdb_id, force)

                # Print result
                pdb_id, success, error = result
                if success:
                    click.echo(f"✓ {pdb_id}: Successfully uploaded")
                else:
                    click.echo(f"✗ {pdb_id}: {error}")

                queue.task_done()

            except asyncio.CancelledError:
                break

    async def main():
        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Get list of structures
        structures = GlobalOps.list_profiles()
        total = len(structures)

        click.echo(f"Found {total} structures to upload")

        # Create queue and workers
        queue = asyncio.Queue()
        worker_tasks = []

        # Start workers
        for _ in range(workers):
            task = asyncio.create_task(upload_worker(queue, adapter, force))
            worker_tasks.append(task)

        # Add structures to queue
        for pdb_id in structures:
            await queue.put(pdb_id)

        # Wait for all uploads to complete
        await queue.join()

        # Cancel workers
        for task in worker_tasks:
            task.cancel()
        await asyncio.gather(*worker_tasks, return_exceptions=True)

    try:
        # Run the async main function
        asyncio.run(main())
        click.echo("\nUpload complete!")

    except Exception as e:
        click.echo(f"Error during upload: {str(e)}", err=True)
        raise click.Abort()


@db.command()
@click.option("--workers", "-w", type=int, default=4, help="Number of parallel workers")
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show which structures would be uploaded without uploading",
)
@click.option("--uri", envvar="NEO4J_URI", help="Neo4j database URI")
@click.option("--user", envvar="NEO4J_USER", help="Neo4j username")
@click.option("--password", envvar="NEO4J_PASSWORD", help="Neo4j password")
@click.option(
    "--database", envvar="NEO4J_DATABASE", help="Neo4j database name", default="neo4j"
)
def upload_missing(
    workers: int,
    dry_run: bool,
    uri: Optional[str],
    user: Optional[str],
    password: Optional[str],
    database: str,
) -> None:
    """
    Upload structures that exist in RCSB but not in the local database.

    Uses GlobalOps.status_vs_rcsb() to identify missing structures and uploads them in parallel.
    """

    async def upload_structure(
        adapter: Neo4jAdapter, pdb_id: str
    ) -> tuple[str, bool, Optional[str]]:
        """Upload a single structure and return result"""
        try:
            adapter.add_total_structure(pdb_id, disable_exists_check=True)
            return pdb_id, True, None
        except Exception as e:
            return pdb_id, False, str(e)

    async def upload_worker(queue: asyncio.Queue, adapter: Neo4jAdapter) -> None:
        """Worker to process structures from the queue"""
        while True:
            try:
                pdb_id = await queue.get()
                result = await upload_structure(adapter, pdb_id)

                # Print result
                pdb_id, success, error = result
                if success:
                    click.echo(f"✓ {pdb_id}: Successfully uploaded")
                else:
                    click.echo(f"✗ {pdb_id}: {error}")

                queue.task_done()

            except asyncio.CancelledError:
                break

    async def main():
        # Get list of missing structures
        db_entries = Neo4jReader(
            Neo4jAdapter(
                uri=NEO4J_URI,
                user=NEO4J_USER,
                password=NEO4J_PASSWORD,
                current_db=NEO4J_CURRENTDB,
            )
        ).all_ids()
        missing_structures = GlobalOps.missing_db_entries(db_entries)
        total = len(missing_structures)

        if total == 0:
            click.echo("No missing structures found!")
            return

        for struct in missing_structures:
            click.echo(f"- {struct}")
        click.echo(
            f"Found {total} structures in RCSB that are not in the local database:"
        )

        if dry_run:
            click.echo("\nDry run - no structures were uploaded")
            return

        if not click.confirm(
            f"\nDo you want to proceed with uploading {total} structures?"
        ):
            click.echo("Upload cancelled")
            return

        # Initialize database adapter
        adapter = Neo4jAdapter(
            uri=NEO4J_URI,
            user=NEO4J_USER,
            password=NEO4J_PASSWORD,
            current_db=NEO4J_CURRENTDB,
        )

        # Create queue and workers
        queue = asyncio.Queue()
        worker_tasks = []

        # Start workers
        for _ in range(workers):
            task = asyncio.create_task(upload_worker(queue, adapter))
            worker_tasks.append(task)

        with click.progressbar(
            missing_structures, label="Queueing structures", length=total
        ) as structures:
            for pdb_id in structures:
                await queue.put(pdb_id)

        # Wait for all uploads to complete
        await queue.join()

        # Cancel workers
        for task in worker_tasks:
            task.cancel()
        await asyncio.gather(*worker_tasks, return_exceptions=True)

    try:
        # Run the async main function
        asyncio.run(main())
        click.echo("\nUpload complete!")

    except Exception as e:
        click.echo(f"Error during upload: {str(e)}", err=True)
        raise click.Abort()


@db.group()
@click.pass_context
def instance(ctx):
    """Manage Neo4j database instances"""
    pass


@instance.command()
@click.option(
    "--name", prompt="Database name", help="Name for the new database instance"
)
@click.option(
    "--initialize/--no-initialize",
    default=True,
    help="Initialize with constraints and basic data",
)
@click.option(
    "--force", "-f", is_flag=True, help="Force creation even if database exists"
)
def create(name: str, initialize: bool, force: bool):
    """Create a new Neo4j database instance.

    Examples:\n
    \b
    # Create a new database with prompts
    ribd.py db instance create

    \b
    # Create a specific database non-interactively
    ribd.py db instance create --name ribosome_test

    \b
    # Create without initialization
    ribd.py db instance create --name ribosome_test --no-initialize
    """
    try:
        # Connect to system database to manage instances
        adapter = Neo4jAdapter(
            NEO4J_URI,
            NEO4J_USER,
            "system",  # Use system database for management
            NEO4J_PASSWORD,
        )

        # Check if database exists
        if not force:
            with adapter.driver.session() as session:
                result = session.run("SHOW DATABASES WHERE name = $name", name=name)
                if result.single():
                    if not click.confirm(
                        f'Database "{name}" already exists. Override?'
                    ):
                        click.echo("Cancelled.")
                        return

        # Create database
        with adapter.driver.session() as session:
            session.run(f"CREATE DATABASE {name} IF NOT EXISTS")
            click.echo(f"Created database: {name}")

        if initialize:
            # Switch adapter to new database
            adapter = Neo4jAdapter(
                NEO4J_URI, NEO4J_USER, name, NEO4J_PASSWORD  # Use new database
            )

            with click.progressbar(
                length=3, label="Initializing database", show_pos=True
            ) as bar:
                # Initialize constraints
                adapter.init_constraints()
                bar.update(1)

                # Initialize polymer classes
                adapter.init_polymer_classes()
                bar.update(1)

                # Initialize phylogenies
                adapter.init_phylogenies()
                bar.update(1)

            click.echo("\nDatabase initialized successfully!")

    except Exception as e:
        click.echo(f"Error creating database: {str(e)}", err=True)
        raise


@instance.command()
def list():
    """List all Neo4j database instances"""
    try:
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, "system", NEO4J_PASSWORD)

        with adapter.driver.session() as session:
            result = session.run("SHOW DATABASES")
            databases = result.data()

            if not databases:
                click.echo("No databases found")
                return

            click.echo("\nAvailable databases:")
            for db in databases:
                status = "🟢" if db["currentStatus"] == "online" else "🔴"
                click.echo(f"{status} {db['name']:<20} ({db['currentStatus']})")

    except Exception as e:
        click.echo(f"Error listing databases: {str(e)}", err=True)


@instance.command(name="delete")
@click.argument("name", required=True)
@click.option("--force", "-f", is_flag=True, help="Force deletion without confirmation")
@click.pass_context
def delete_instance(ctx, name, force):
    """Delete a Neo4j database instance and all its data.

    Examples:\n
    \b
    # Delete an instance with confirmation
    ribd.py db instance delete ribosome2

    \b
    # Force delete without confirmation
    ribd.py db instance delete --force ribosome2
    """
    try:
        # Connect to system database for instance management
        adapter = Neo4jAdapter(NEO4J_URI, NEO4J_USER, "system", NEO4J_PASSWORD)

        with adapter.driver.session() as session:
            # Check if instance exists
            result = session.run("SHOW DATABASES")
            databases = [record["name"] for record in result]

            if name not in databases:
                click.echo(f"Database instance '{name}' does not exist.", err=True)
                return

            # Check if it's the default database
            if name == "neo4j":
                click.echo("Cannot delete the default 'neo4j' database.", err=True)
                return

            # Get instance status
            result = session.run("SHOW DATABASE $name", name=name)
            status = result.single()["currentStatus"]

            if not force:
                click.echo(f"\nDatabase: {name}")
                click.echo(f"Status: {status}")
                click.confirm(
                    "Are you sure you want to delete this database and ALL its data?",
                    abort=True,
                )

            # Stop the database if it's running
            if status != "offline":
                click.echo("Stopping database...")
                session.run("STOP DATABASE $name", name=name)

            # Drop the database
            click.echo("Deleting database...")
            session.run("DROP DATABASE $name IF EXISTS", name=name)

            # Verify deletion
            result = session.run("SHOW DATABASES")
            databases = [record["name"] for record in result]

            if name not in databases:
                click.echo(f"Successfully deleted database instance '{name}'")
            else:
                click.echo(f"Failed to delete database instance '{name}'", err=True)

    except Exception as e:
        click.echo(f"Failed to delete database instance: {str(e)}", err=True)


@etl.command()
def bootstrap():
    """Seed the local NCBI taxonomy database and other necessary binary resources."""
    from ribctl.lib.libtax import get_ncbi
    click.echo("Checking NCBI taxonomy database...")
    try:
        ncbi = get_ncbi()
        # This force-checks if the DB is actually queryable
        ncbi.get_taxid_translator([9606])
        click.echo("✓ NCBI taxonomy database is ready.")
    except Exception as e:
        click.echo("NCBI database missing or corrupt. Initializing...")
        from ete3 import NCBITaxa
        ncbi = NCBITaxa(dbfile=os.environ.get("NCBI_TAXA_SQLITE"))
        ncbi.update_taxonomy_database()
        click.echo("✓ NCBI taxonomy database initialized.")

if __name__ == "__main__":
    cli(obj={})


```

ribctl/ribosome_ops.py
```py
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.asset_manager.assets_structure import StructureAssets
from Bio.PDB.Structure import Structure
from ribctl.lib.types.polymer.base import CytosolicRNAClass, MitochondrialRNAClass
from ribctl.lib.schema.types_ribosome import ( RNA, PTCInfo, Polymer, PolymerClass,  RibosomeStructure, RibosomeStructureMetadata, )
from ribctl import RIBETL_DATA

class RibosomeOps:
    rcsb_id: str
    assets: StructureAssets

    def __init__(self, rcsb_id: str) -> None:
        if not RIBETL_DATA:
            raise Exception("RIBETL_DATA environment variable not set. Cannot access assets." )
        self.rcsb_id = rcsb_id.upper()
        self.assets  = StructureAssets(self.rcsb_id)

    @property
    def taxid(self)->int:
        return self.assets.profile().src_organism_ids[0]

    @property
    def profile(self) -> RibosomeStructure:
        return self.assets.profile()


    def first_assembly_auth_asym_ids(self) -> list[str]:
        """
        Returns a list of auth_asym_ids belonging to the first assembly in the structure.
        If no assembly map is present, raises an exception.
        
        Returns:
            list[str]: List of auth_asym_ids from the first assembly
        
        Raises:
            Exception: If no assembly map is found in the structure
        """
        if not self.assets.profile().assembly_map:
            raise Exception("No assembly map found in structure")
        
        # Get the first assembly from the assembly map
        if not self.assets.profile().assembly_map:
            raise Exception("No assembly map found in structure")
        first_assembly = self.assets.profile().assembly_map[0]
        auth_asym_ids = []
        
        # Get auth_asym_ids from polymer instances
        if first_assembly.polymer_entity_instances:
            auth_asym_ids.extend([
                instance.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id 
                for instance in first_assembly.polymer_entity_instances
            ])
        
        # Get auth_asym_ids from nonpolymer instances if they exist
        if first_assembly.nonpolymer_entity_instances:
            auth_asym_ids.extend([
                instance.rcsb_nonpolymer_entity_instance_container_identifiers.auth_asym_id
                for instance in first_assembly.nonpolymer_entity_instances
            ])
        
        return auth_asym_ids

    def get_biopython_chain_by_polymer_class(self, polymer_class: PolymerClass) -> Chain:
        model = self.assets.biopython_structure()[0]
        poly = self.get_poly_by_polyclass(polymer_class)
        if poly == None:
            raise KeyError("No polymer found with class: {}".format(polymer_class))
        return model.child_dict[poly.auth_asym_id]

    def get_biopython_chain_by_auth_asym_id(self, auth_asym_id: str) -> Chain:
        model = self.assets.biopython_structure()[0]
        return model.child_dict[auth_asym_id]

    def nomenclature_table(self, verbose: bool = False) -> dict[str, dict]:
        prof = self.assets.profile()
        m    = {}

        for p in prof.other_polymers:
            m[p.auth_asym_id] = {
                "nomenclature": list(map(lambda x: x.value, p.nomenclature)),
            }

            if verbose:
                m[p.auth_asym_id].update(
                    {
                        "entity_poly_strand_id": p.entity_poly_strand_id,
                        "rcsb_pdbx_description": p.rcsb_pdbx_description,
                    }
                )
        for prot in prof.proteins:
            m[prot.auth_asym_id] = {
                "nomenclature": list(map(lambda x: x.value, prot.nomenclature)),
            }
            if verbose:
                m[prot.auth_asym_id].update(
                    {
                        "entity_poly_strand_id": prot.entity_poly_strand_id,
                        "rcsb_pdbx_description": prot.rcsb_pdbx_description,
                    }
                )
        if prof.rnas != None:
            for rna in prof.rnas:
                m[rna.auth_asym_id] = {
                    "nomenclature": list(map(lambda x: x.value, rna.nomenclature)),
                }
                if verbose:
                    m[rna.auth_asym_id].update(
                        {
                            "entity_poly_strand_id": rna.entity_poly_strand_id,
                            "rcsb_pdbx_description": rna.rcsb_pdbx_description,
                        }
                    )

        return m

    def get_taxids(self) -> tuple[list[int], list[int]]:
        p = self.assets.profile()
        return (p.src_organism_ids, p.host_organism_ids)

    def get_poly_by_auth_asym_id( self, auth_asym_id: str ) -> Polymer :
        profile = self.assets.profile()
        for chain in [ *profile.proteins, *profile.rnas, *profile.other_polymers]:
            if chain.auth_asym_id == auth_asym_id:
                return chain
        raise KeyError("No chain found with auth_asym_id: {}".format(auth_asym_id))

    def get_poly_by_polyclass( self, class_: PolymerClass, assembly: int = 0 ) ->Polymer | None:
        """@assembly here stands to specify which of the two or more models the rna comes from
        in the case that a structure contains multiple models (ex. 4V4Q XRAY)"""

        profile = self.assets.profile()
        polymer:Polymer
        for polymer in [*profile.rnas, *profile.other_polymers, *profile.proteins]: 
            if class_ in  polymer.nomenclature and polymer.assembly_id == assembly and polymer.entity_poly_seq_length > 30:
                return polymer

    def get_LSU_rRNA(self, assembly: int = 0) -> RNA:
        """retrieve the largest rRNA sequence in the structure
        @returns (seq, auth_asym_id, rna_type)
        """

        rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_23S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_25S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(CytosolicRNAClass.rRNA_28S, assembly)
        if rna == None:
            rna = self.get_poly_by_polyclass(MitochondrialRNAClass.mtrRNA16S, assembly)
        if rna == None:
            raise Exception("No LSU rRNA found in structure")
        else:
            return rna

    @staticmethod
    def biopython_chain_get_seq(
        struct: Structure,
        auth_asym_id: str,
        protein_rna: typing.Literal["protein", "rna"],
        sanitized: bool = False,
    ) -> str:

        chain3d = struct.child_dict[0].child_dict[auth_asym_id]
        ress    = chain3d.child_list

        seq = ""
        for i in ress:
            if i.resname not in AMINO_ACIDS_3_TO_1_CODE.keys():
                print("Unknown residue:", i.resname)
                continue

            if protein_rna == "rna":
                seq += i.resname
            else:
                seq += AMINO_ACIDS_3_TO_1_CODE[i.resname]

        return seq
```

ribctl/lib/types/polymer/types.py
```py
from typing import Union
from ribctl.lib.types.polymer.base import (
    CytosolicProteinClass,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    MitochondrialProteinClass,
    MitochondrialRNAClass,
    tRNA,
)

ProteinClass         = Union[MitochondrialProteinClass, CytosolicProteinClass]
LifecycleFactorClass = Union[ElongationFactorClass, InitiationFactorClass]
PolypeptideClass     = Union[LifecycleFactorClass, ProteinClass]
PolynucleotideClass  = Union[CytosolicRNAClass, MitochondrialRNAClass, tRNA]
PolymerClass         = Union[PolynucleotideClass, PolypeptideClass]

```

ribctl/lib/types/polymer/hierarchies.py
```py
from enum import Enum
from typing import  Union, Type, Iterator, TypeVar, Generic

from ribctl.lib.types.polymer.base import CytosolicProteinClass, CytosolicRNAClass, ElongationFactorClass, InitiationFactorClass, MitochondrialProteinClass, MitochondrialRNAClass, tRNA

E = TypeVar("E", bound=Enum)


class PolymerHierarchy(Generic[E]):
    """
    Manages a hierarchical collection of related Enum classes.
    Provides type-safe iteration, containment checking, and classification.
    """

    def __init__(self, *enum_classes: Type[E], name: None | str = None):
        self.enum_classes = enum_classes
        self.name = name or "+".join(cls.__name__ for cls in enum_classes)

        # Pre-compute all valid values for faster lookup
        self._all_values = {
            member.value: (cls, member) for cls in enum_classes for member in cls
        }

    def __iter__(self) -> Iterator[E]:
        """Allows iteration over all members across all enum classes."""
        return (member for cls in self.enum_classes for member in cls)

    def __contains__(self, item: Union[str, E]) -> bool:
        """Enables 'in' operator for both string values and enum members."""
        if isinstance(item, str):
            return item in self._all_values
        return item in self._all_values.values()

    def get_type(self, value: str) -> Type[E] | None:
        """Returns the enum class type for a given value."""

        if value in self._all_values:
            return self._all_values[value][0]
        return None

    def get_member(self, value: str) -> E | None:
        """Returns the enum member for a given value."""
        if value in self._all_values:
            return self._all_values[value][1]
        return None

    def is_of_type(self, value: Union[str, E], enum_type: Type[E]) -> bool:
        """Checks if a value belongs to a specific enum type."""
        if isinstance(value, str):
            cls = self.get_type(value)
            return cls == enum_type if cls else False
        return isinstance(value, enum_type)


Proteins         = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, name="Proteins" )
LifecycleFactors = PolymerHierarchy( ElongationFactorClass, InitiationFactorClass, name="LifecycleFactors" )
Polypeptides     = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, name="Polypeptides", )
Polynucleotides  = PolymerHierarchy( CytosolicRNAClass, MitochondrialRNAClass, tRNA, name="Polynucleotides" )
Polymers         = PolymerHierarchy( MitochondrialProteinClass, CytosolicProteinClass, ElongationFactorClass, InitiationFactorClass, CytosolicRNAClass, MitochondrialRNAClass, tRNA, name="Polymers" )
```

ribctl/lib/types/polymer/base.py
```py
from enum import Enum

class PolymerEnumBase(str, Enum):
    def __repr__(self):
        return self.value

class tRNA(PolymerEnumBase):
    tRNA = "tRNA"

class MitochondrialProteinClass(PolymerEnumBase):
    # mSSU
    bS1m  = "bS1m"
    uS2m  = "uS2m"
    uS3m  = "uS3m"
    uS4m  = "uS4m"
    uS5m  = "uS5m"
    bS6m  = "bS6m"
    uS7m  = "uS7m"
    uS8m  = "uS8m"
    uS9m  = "uS9m"
    uS10m = "uS10m"
    uS11m = "uS11m"
    uS12m = "uS12m"
    uS13m = "uS13m"
    uS14m = "uS14m"
    uS15m = "uS15m"
    bS16m = "bS16m"
    uS17m = "uS17m"
    bS18m = "bS18m"
    uS19m = "uS19m"
    bS21m = "bS21m"
    mS22  = "mS22"
    mS23  = "mS23"
    mS25  = "mS25"
    mS26  = "mS26"
    mS27  = "mS27"
    mS29  = "mS29"
    mS31  = "mS31"
    mS33  = "mS33"
    mS34  = "mS34"
    mS35  = "mS35"
    mS37  = "mS37"
    mS38  = "mS38"
    mS39  = "mS39"
    mS40  = "mS40"
    mS41  = "mS41"
    mS42  = "mS42"
    mS43  = "mS43"
    mS44  = "mS44"
    mS45  = "mS45"
    mS46  = "mS46"
    mS47  = "mS47"

    # mLSU
    uL1m  = "uL1m"
    uL2m  = "uL2m"
    uL3m  = "uL3m"
    uL4m  = "uL4m"
    uL5m  = "uL5m"
    uL6m  = "uL6m"
    bL9m  = "bL9m"
    uL10m = "uL10m"
    uL11m = "uL11m"
    bL12m = "bL12m"
    uL13m = "uL13m"
    uL14m = "uL14m"
    uL15m = "uL15m"
    uL16m = "uL16m"
    bL17m = "bL17m"
    uL18m = "uL18m"
    bL19m = "bL19m"
    bL20m = "bL20m"
    bL21m = "bL21m"
    uL22m = "uL22m"
    uL23m = "uL23m"
    uL24m = "uL24m"
    bL27m = "bL27m"
    bL28m = "bL28m"
    uL29m = "uL29m"
    uL30m = "uL30m"
    bL31m = "bL31m"
    bL32m = "bL32m"
    bL33m = "bL33m"
    bL34m = "bL34m"
    bL35m = "bL35m"
    bL36m = "bL36m"
    mL37  = "mL37"
    mL38  = "mL38"
    mL39  = "mL39"
    mL40  = "mL40"
    mL41  = "mL41"
    mL42  = "mL42"
    mL43  = "mL43"
    mL44  = "mL44"
    mL45  = "mL45"
    mL46  = "mL46"
    mL48  = "mL48"
    mL49  = "mL49"
    mL50  = "mL50"
    mL51  = "mL51"
    mL52  = "mL52"
    mL53  = "mL53"
    mL54  = "mL54"
    mL57  = "mL57"
    mL58  = "mL58"
    mL59  = "mL59"
    mL60  = "mL60"
    mL61  = "mL61"
    mL62  = "mL62"
    mL63  = "mL63"
    mL64  = "mL64"
    mL65  = "mL65"
    mL66  = "mL66"
    mL67  = "mL67"

class CytosolicProteinClass(PolymerEnumBase):
    # SSU

    
    bS1   = "bS1"
    eS1   = "eS1"
    uS2   = "uS2"
    uS3   = "uS3"
    uS4   = "uS4"
    eS4   = "eS4"
    uS5   = "uS5"
    bS6   = "bS6"
    eS6   = "eS6"
    uS7   = "uS7"
    eS7   = "eS7"
    uS8   = "uS8"
    eS8   = "eS8"
    uS9   = "uS9"
    uS10  = "uS10"
    eS10  = "eS10"
    uS11  = "uS11"
    uS12  = "uS12"
    eS12  = "eS12"
    uS13  = "uS13"
    uS14  = "uS14"
    uS15  = "uS15"
    bS16  = "bS16"
    uS17  = "uS17"
    eS17  = "eS17"
    bS18  = "bS18"
    uS19  = "uS19"
    eS19  = "eS19"
    bS20  = "bS20"
    bS21  = "bS21"
    bTHX  = "bTHX"
    eS21  = "eS21"
    eS24  = "eS24"
    eS25  = "eS25"
    eS26  = "eS26"
    eS27  = "eS27"
    eS28  = "eS28"
    eS30  = "eS30"
    eS31  = "eS31"
    RACK1 = "RACK1"
    # LSU
    uL1  = "uL1"
    uL2  = "uL2"
    uL3  = "uL3"
    uL4  = "uL4"
    uL5  = "uL5"
    uL6  = "uL6"
    eL6  = "eL6"
    eL8  = "eL8"
    bL9  = "bL9"
    uL10 = "uL10"
    uL11 = "uL11"
    bL12 = "bL12"
    uL13 = "uL13"
    eL13 = "eL13"
    uL14 = "uL14"
    eL14 = "eL14"
    uL15 = "uL15"
    eL15 = "eL15"
    uL16 = "uL16"
    bL17 = "bL17"
    uL18 = "uL18"
    eL18 = "eL18"
    bL19 = "bL19"
    eL19 = "eL19"
    bL20 = "bL20"
    eL20 = "eL20"
    bL21 = "bL21"
    eL21 = "eL21"
    uL22 = "uL22"
    eL22 = "eL22"
    uL23 = "uL23"
    uL24 = "uL24"
    eL24 = "eL24"
    bL25 = "bL25"
    bL27 = "bL27"
    eL27 = "eL27"
    bL28 = "bL28"
    eL28 = "eL28"
    uL29 = "uL29"
    eL29 = "eL29"
    uL30 = "uL30"
    eL30 = "eL30"
    bL31 = "bL31"
    eL31 = "eL31"
    bL32 = "bL32"
    eL32 = "eL32"
    bL33 = "bL33"
    eL33 = "eL33"
    bL34 = "bL34"
    eL34 = "eL34"
    bL35 = "bL35"
    bL36 = "bL36"
    eL36 = "eL36"
    eL37 = "eL37"
    eL38 = "eL38"
    eL39 = "eL39"
    eL40 = "eL40"
    eL41 = "eL41"
    eL42 = "eL42"
    eL43 = "eL43"
    P1P2 = "P1P2"

class MitochondrialRNAClass(PolymerEnumBase):
    mtrRNA12S = "mt12SrRNA"  # mitochondrial
    mtrRNA16S = "mt16SrRNA"  # mitochondrial
    
class CytosolicRNAClass(PolymerEnumBase):
    rRNA_5S   = "5SrRNA"  #  bacterial or eykaryotic
    rRNA_16S  = "16SrRNA"  #  c-bacterial or mitochondrial
    rRNA_23S  = "23SrRNA"  # bacterial
    rRNA_25S  = "25SrRNA"  # plants
    rRNA_5_8S = "5.8SrRNA"  # eukaryotic
    rRNA_18S  = "18SrRNA"  # eukaryotic
    rRNA_28S  = "28SrRNA"  # eukaryotic

class ElongationFactorClass(PolymerEnumBase):
    # Eukaryotic
    eEF1A = "eEF1A"
    eEF1B = "eEF1B"
    eFSec = "eFSec"
    eEF2  = "eEF2"
    mtEF4 = "mtEF4"
    eIF5A = "eIF5A"
    eEF3  = "eEF3"
    # Bacterial
    EF_Tu = "EF-Tu"
    EF_Ts = "EF-Ts"
    SelB  = "SelB"
    EF_G  = "EF-G"
    EF4   = "EF4"
    EF_P  = "EF-P"
    Tet_O = "Tet_O"
    Tet_M = "Tet_M"
    RelA  = "RelA"
    BipA  = "BipA"
    # Archaeal
    aEF1A = "aEF1A"
    aEF2  = "aEF2"

class InitiationFactorClass(PolymerEnumBase):
    #!Eukaryotic
    eIF1  = "eIF1"
    eIF1A = "eIF1A"

    eIF2_alpha = "eIF2_alpha"
    eIF2_beta  = "eIF2_beta"
    eIF2_gamma = "eIF2_gamma"

    eIF2B_alpha   = "eIF2B_alpha"
    eIF2B_beta    = "eIF2B_beta"
    eIF2B_gamma   = "eIF2B_gamma"
    eIF2B_delta   = "eIF2B_delta"
    eIF2B_epsilon = "eIF2B_epsilon"

    eIF3_subunitA = "eIF3_subunitA"
    eIF3_subunitB = "eIF3_subunitB"
    eIF3_subunitC = "eIF3_subunitC"
    eIF3_subunitD = "eIF3_subunitD"
    eIF3_subunitE = "eIF3_subunitE"
    eIF3_subunitF = "eIF3_subunitF"
    eIF3_subunitG = "eIF3_subunitG"
    eIF3_subunitH = "eIF3_subunitH"
    eIF3_subunitI = "eIF3_subunitI"
    eIF3_subunitJ = "eIF3_subunitJ"
    eIF3_subunitK = "eIF3_subunitK"
    eIF3_subunitL = "eIF3_subunitL"
    eIF3_subunitM = "eIF3_subunitM"

    eIF4F_4A = "eIF4F_4A"
    eIF4F_4G = "eIF4F_4G"
    eIF4F_4E = "eIF4F_4E"

    eIF4B = "eIF4B"
    eIF5B = "eIF5B"
    eIF5  = "eIF5"

    #!Bacterial
    IF1 = "IF1"
    IF2 = "IF2"
    IF3 = "IF3"

    #!Archaeal
    aIF_1A       = "aIF1A"
    aIF_2_alpha  = "aIF2_alpha"
    aIF_2_beta   = "aIF2_beta"
    aIF_2_gamma  = "aIF2_gamma"
    aIF_2B_alpha = "aIF2B_alpha"
    aIF_2B_beta  = "aIF2B_beta"
    aIF_2B_delta = "aIF2B_delta"
    aIF5A        = "aIF5A"
    aIF5B        = "aIF5B"

```

ribctl/lib/schema/types_ribosome.py
```py
from typing import Dict, Optional
from enum import Enum
import typing
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.PDB.Residue import Residue
from pydantic import BaseModel
from typing import Optional
from ribctl.lib.types.polymer import PolymerClass
from ribctl.lib.schema.primitives import AMINO_ACIDS, NUCLEOTIDES


class Polymer(BaseModel):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    def to_SeqRecord(self) -> SeqRecord:
        return SeqRecord(
            seq         = Seq(self.entity_poly_seq_one_letter_code_can),
            id          = f"{self.src_organism_ids[0]}",
            description = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id),
            name        = '{}.{}'.format(self.parent_rcsb_id,self.auth_asym_id)
        )

    assembly_id: int

    asym_ids    : list[str]
    auth_asym_id: str

    parent_rcsb_id: str

    src_organism_names : list[str]
    host_organism_names: list[str]

    src_organism_ids   : list[int]
    host_organism_ids  : list[int]

    rcsb_pdbx_description              : Optional[str] = None

    entity_poly_strand_id              : str
    entity_poly_seq_one_letter_code    : str
    entity_poly_seq_one_letter_code_can: str
    entity_poly_seq_length             : int
    entity_poly_polymer_type           : str
    entity_poly_entity_type            : str

    nomenclature                       : list[PolymerClass]

class Protein(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    @staticmethod
    def from_polymer(p: Polymer, **kwargs):
        if kwargs["pfams"] != None and len(kwargs["pfams"]) > 0:

            pfam_comments     = list( set([pfam["rcsb_pfam_comment"] for pfam in kwargs["pfams"]]) )
            pfam_descriptions = list( set([pfam["rcsb_pfam_description"] for pfam in kwargs["pfams"]]) )
            pfam_accessions   = list( set([pfam["rcsb_pfam_accession"] for pfam in kwargs["pfams"]]) )

        else:
            pfam_comments = []
            pfam_descriptions = []
            pfam_accessions = []

        return Protein(
            **{
                **p.model_dump(),
                "pfam_accessions": pfam_accessions,
                "pfam_comments": pfam_comments,
                "pfam_descriptions": pfam_descriptions,
                "uniprot_accession": [entry["rcsb_id"] for entry in kwargs["uniprots"]]
                if kwargs["uniprots"] != None and len(kwargs["uniprots"]) > 0
                else [],
            }
        )

    pfam_accessions  : list[str]
    pfam_comments    : list[str]
    pfam_descriptions: list[str]
    uniprot_accession: list[str]

    def to_polymer(self) -> Polymer:
        return Polymer(**self.model_dump())

class RNA(Polymer):
    def __hash__(self):
        return hash(self.auth_asym_id + self.parent_rcsb_id)

    # pass

class NonpolymericLigand(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    # def metadatum(self) -> NonpolymericLigandMetadatum:
    #     return NonpolymericLigandMetadatum(**self.model_dump())

    class NonpolymerComp(BaseModel):

        class Drugbank(BaseModel):
            class DrugbankInfo(BaseModel):
                cas_number: Optional[str] =None
                description: Optional[str] =None

            class DrugbankContainerIdentifiers(BaseModel):
                drugbank_id: str

            drugbank_container_identifiers: DrugbankContainerIdentifiers
            drugbank_info: DrugbankInfo

        class RcsbChemCompTarget(BaseModel):
            interaction_type: Optional[str] = None
            name: Optional[str] =None
            provenance_source: Optional[str] =None
            reference_database_accession_code: Optional[str] =None
            reference_database_name: Optional[str] =None

        drugbank             : Optional[Drugbank]=None
        rcsb_chem_comp_target: Optional[list[RcsbChemCompTarget]]=None

    chemicalId    : str
    chemicalName  : str
    formula_weight: Optional[float] =None

    pdbx_description   : str
    number_of_instances: int

    nonpolymer_comp: Optional[NonpolymerComp] = None

    SMILES       : Optional[str] =None
    SMILES_stereo: Optional[str] =None
    InChI        : Optional[str] =None
    InChIKey     : Optional[str] =None

class ResidueSummary(BaseModel): 
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }

    label_seq_id : typing.Optional[int] = None
    label_comp_id: typing.Optional[str] = None
    auth_asym_id : str
    auth_seq_id  : int
    rcsb_id      : str
    full_id      : typing.Optional[tuple[str, int, str, tuple[str, int, str]]]


    @staticmethod
    def is_canonical(resname:str):
        return resname in [*AMINO_ACIDS.keys(), *NUCLEOTIDES]

    @staticmethod
    def three_letter_code_to_one(resname: str):
        if resname in AMINO_ACIDS:
            return AMINO_ACIDS[resname]["one_letter_code"]
        elif resname in NUCLEOTIDES:
                return resname
        else:
            return '-'

    @staticmethod
    def one_letter_code_to_three(resname: str):
        if resname in [*map(lambda x: x[1]['one_letter_code'], AMINO_ACIDS.items())]:
            for tlk, d in AMINO_ACIDS.items():
                if d["one_letter_code"] == resname:
                    return tlk
        elif resname in NUCLEOTIDES:
                return resname
        else:
            return '-'

    def __hash__(self):
        return hash( self.get_resname() if self.get_resname() is not None else "" + str(self.get_seqid()) + self.get_parent_auth_asym_id() )

    def get_resname(self):
        return self.label_comp_id

    def get_seqid(self):
        (structure_id, model_id, chain_id, _) = self.full_id
        (hetero, seqid, insertion_code)       = _
        return seqid

    def get_parent_auth_asym_id(self):
        (structure_id, assembly_id, chain_id, _) = self.full_id
        return chain_id

    @staticmethod
    def from_biopython_residue(r: Residue):

        (structure_id, model_id, chain_id, _) = r.get_full_id()
        (hetero, seqid, insertion_code) = _

        return ResidueSummary(
            auth_seq_id   = seqid,
            label_seq_id  = None,
            label_comp_id = r.get_resname(),
            auth_asym_id  = chain_id,
            full_id       = r.get_full_id(),
            rcsb_id       = r.get_full_id()[0]
        )

class AssemblyInstancesMap(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    """
    This basically specifies which assembly an instnace of a polymer or a nonpolymer belongs to.
    Certain PDB structures come with more than a single physical model/assembly packaged in the file,
    hence every chain and many ligands might be present in 2 or more instances.

    The RNA/Protein/Ligand suffices to characterizes all instances, yet, to resolve duplicate chains
    precisely in space, this information is needed.

    assemblies{
    rcsb_id
        nonpolymer_entity_instances{
          rcsb_nonpolymer_entity_instance_container_identifiers{
            entity_id
            comp_id
            auth_asym_id
            rcsb_id
            auth_seq_id
          }
        }
        polymer_entity_instances{
           rcsb_polymer_entity_instance_container_identifiers {
            entity_id
            asym_id
            auth_asym_id
            entry_id
            entity_id
          }
        }
      }
    """

    class NonpolymerEntityInstance(BaseModel):
        class NonpolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str
            auth_seq_id: str

        rcsb_nonpolymer_entity_instance_container_identifiers: NonpolymerEntityInstanceContainerIdentifiers

    class PolymerEntityInstance(BaseModel):
        class PolymerEntityInstanceContainerIdentifiers(BaseModel):
            entity_id: str
            auth_asym_id: str
        rcsb_polymer_entity_instance_container_identifiers: PolymerEntityInstanceContainerIdentifiers

    rcsb_id                    : str  # ex. 5AFI-1
    nonpolymer_entity_instances: Optional[list[NonpolymerEntityInstance]] =None
    polymer_entity_instances   : list[PolymerEntityInstance]

class NomenclatureItem(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    nomenclature: list[str]

class NomenclatureTable(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    __pydantic_root_model__: Dict[str, NomenclatureItem]

class PTCInfo(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    location: list[float]
    residues: list[ResidueSummary]

class ConstrictionSite(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }
    location: list[float]

class RibosomeStructureMetadata(BaseModel):
    model_config = {
        "json_encoders": {
            Enum: lambda v: v.value
        }
    }

    def get_polymers_by_assembly(self)->dict[str,list[str]]:
        if self.assembly_map:
            _ = {}
            for assembly in self.assembly_map:
                _[assembly.rcsb_id] = [poly.rcsb_polymer_entity_instance_container_identifiers.auth_asym_id for poly in assembly.polymer_entity_instances]
            return _
        else:
            raise Exception("No assembly map found")

    def get_nomenclature_map(self) -> dict[str, list[PolymerClass]]:
        "Return a map from auth_asym_id to nomenclatures of all polymers in the structure"
        _ = {}

        for prp in self.proteins:
            _[prp.auth_asym_id] = [ x.name for x in  prp.nomenclature]

        for prna in self.rnas:
            _[prna.auth_asym_id] =[ y.name for y in  prna.nomenclature]

        return _




    def __hash__(self):
        return hash(self.rcsb_id)

    @staticmethod
    def model_validate_partial(data: dict, polymer_fields: bool = False):
        if polymer_fields:
            return RibosomeStructure.model_validate(data)
        else:
            # Create a copy of the data without the detailed fields
            filtered_data = {k: v for k, v in data.items() if k not in ['proteins', 'rnas', 'other_polymers', 'nonpolymeric_ligands']}
            return RibosomeStructure.model_validate(filtered_data)

    rcsb_id   : str
    expMethod : str
    resolution: float

    deposition_date:Optional[ str ] = None

    pdbx_keywords      : Optional[str] = None
    pdbx_keywords_text: Optional[str]  = None

    rcsb_external_ref_id  : list[str]
    rcsb_external_ref_type: list[str]
    rcsb_external_ref_link: list[str]

    citation_year         : None| Optional[int|str] = None
    citation_rcsb_authors: None|Optional[list[str]] = None
    citation_title        :None| Optional[str]      = None
    citation_pdbx_doi     :None| Optional[str]      = None

    src_organism_ids  : list[int]
    src_organism_names: list[str]

    host_organism_ids  : list[int]
    host_organism_names: list[str]

    assembly_map    : Optional[list[AssemblyInstancesMap]]  = None
    mitochondrial   : bool
    subunit_presence: Optional[list[typing.Literal['ssu','lsu']]] = None

class RibosomeStructure(RibosomeStructureMetadata):

    rnas                : list[RNA]
    other_polymers      : list[Polymer]
    nonpolymeric_ligands: list[NonpolymericLigand]
    proteins            : list[Protein]


```

ribctl/lib/schema/primitives.py
```py

AMINO_ACIDS = {
    "ALA": {"one_letter_code": "A", "charge": 0},
    "ARG": {"one_letter_code": "R", "charge": 1},
    "ASN": {"one_letter_code": "N", "charge": 0},
    "ASP": {"one_letter_code": "D", "charge": -1},
    "CYS": {"one_letter_code": "C", "charge": 0},
    "GLU": {"one_letter_code": "E", "charge": -1},
    "GLN": {"one_letter_code": "Q", "charge": 0},
    "GLY": {"one_letter_code": "G", "charge": 0},
    "HIS": {"one_letter_code": "H", "charge": 0},
    "ILE": {"one_letter_code": "I", "charge": 0},
    "LEU": {"one_letter_code": "L", "charge": 0},
    "LYS": {"one_letter_code": "K", "charge": 1},
    "MET": {"one_letter_code": "M", "charge": 0},
    "PHE": {"one_letter_code": "F", "charge": 0},
    "PRO": {"one_letter_code": "P", "charge": 0},
    "SER": {"one_letter_code": "S", "charge": 0},
    "THR": {"one_letter_code": "T", "charge": 0},
    "TRP": {"one_letter_code": "W", "charge": 0},
    "TYR": {"one_letter_code": "Y", "charge": 0},
    "VAL": {"one_letter_code": "V", "charge": 0},
}
NUCLEOTIDES = ["A", "T", "C", "G", "U"]
```


