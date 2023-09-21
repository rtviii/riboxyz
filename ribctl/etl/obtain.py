
import asyncio
from concurrent import futures
import json
import logging
import typing
from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from typing import Optional
from api.logs.loggers import get_updates_logger
from ribctl import AMINO_ACIDS_3_TO_1_CODE
from ribctl.etl.ribosome_assets import Assetlist, RibosomeAssets
from ribctl.lib.tunnel import ptc_resdiues_get, ptc_residues_calculate_midpoint
from ribctl.lib.ribosome_types.types_binding_site import BindingSite
from ribctl.lib.ribosome_types.types_poly_nonpoly_ligand import PolymericFactorClass, RNAClass
from ribctl.lib.mod_extract_bsites import bsite_nonpolymeric_ligand, struct_ligand_ids, struct_polymeric_factor_ids, bsite_polymeric_factor, bsite_polymeric_factor
from ribctl.lib.mod_split_rename import split_rename
from ribctl.etl.etl_pipeline import current_rcsb_structs, ReannotationPipeline, rcsb_single_structure_graphql, query_rcsb_api
from ribctl.lib.mod_render_thumbnail import render_thumbnail
from ribctl.lib.ribosome_types.types_ribosome import RNA, PolymerClass, PolymericFactor, Protein, ProteinClass, RibosomeStructure
from ribctl import RIBETL_DATA
from concurrent.futures import ALL_COMPLETED, Future, ProcessPoolExecutor, ThreadPoolExecutor, wait
logging.getLogger('asyncio').setLevel(logging.WARNING)

logging.basicConfig(
    level=logging.DEBUG,  
    format='%(asctime)s [%(levelname)s] [%(name)s] %(message)s',  
    handlers=[
        logging.StreamHandler(),  
        logging.FileHandler('etl.log')  
    ]
)


import os
async def obtain_assets(rcsb_id: str, assetlist: Assetlist, overwrite: bool = False):
    """Obtain all assets for a given RCSB ID"""

    rcsb_id = rcsb_id.upper()
    # Obtaining assets for 
    logging.debug("※\t\t\t{}".format(rcsb_id))
    assets = RibosomeAssets(rcsb_id)
    assets._verify_dir_exists()

    coroutines = []

    if assetlist.profile:
        coroutines.append(assets._verify_json_profile(overwrite))

    if assetlist.cif:
        coroutines.append(assets._verify_cif(overwrite))

    if assetlist.factors_and_ligands:
        coroutines.append(assets._verify_ligads_and_ligandlike_polys(overwrite))

    if assetlist.ptc_coords:
            ress, auth_asym_id = ptc_resdiues_get(assets.biopython_structure(),  assets.profile().rnas)
            midpoint_coords = ptc_residues_calculate_midpoint(ress, auth_asym_id)
            # residue_labels = [(res.get_resname(), res.id[1]) for res in ress]
            writeout = {
                "site_9_residues"      : [(res.get_resname(), res.get_segid()) for res in ress],
                "LSU_rRNA_auth_asym_id": auth_asym_id,
                "midpoint_coordinates" : midpoint_coords,
                'nomenclature_table'   : assets.nomenclature_table()

            }

            asset_ptc_coords_path = os.path.join(assets._dir_path(),f'{assets.rcsb_id}_PTC_COORDINATES.json')

            with open(asset_ptc_coords_path, 'w') as f:
                logging.debug(f'Writing PTC coordinates to {asset_ptc_coords_path}')
                json.dump(writeout, f)

    if assetlist.cif_modified_and_chains:
        coroutines.append(assets._verify_cif_modified_and_chains(overwrite))

    await asyncio.gather(*coroutines)

def obtain_assets_threadpool(targets: list[str], assetlist: Assetlist, workers: int = 5, get_all: bool = False, overwrite=False):
    """Get all ribosome profiles from RCSB via a threadpool"""
    logger = get_updates_logger()
    if get_all:
        unsynced = sorted(current_rcsb_structs())
    else:
        unsynced = list(map(lambda _: _.upper(), targets))

    tasks: list[Future] = []
    results=[]
    logger.info("Begun downloading ribosome profiles via RCSB")

    def log_commit_result(rcsb_id: str):
        def _(f: Future):
            if not None == f.exception():
                logger.error(rcsb_id + ":" + f.exception().__str__())
            else:
                logger.info(rcsb_id + ": processed successfully.")
        return _

    with ThreadPoolExecutor(max_workers=workers) as executor:
        for rcsb_id in unsynced:
            fut = executor.submit(asyncio.run, obtain_assets(rcsb_id, assetlist, overwrite))
            fut.add_done_callback(log_commit_result(rcsb_id))
            tasks.append(fut)
        for future in futures.as_completed(tasks):
            results.append(future.result())

    logger.info("Finished syncing with RCSB")
