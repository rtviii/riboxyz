import asyncio
import os
from ribctl import RIBETL_DATA
from ribctl.ribosome_assets import Assetlist, RibosomeAssets, obtain_assets, obtain_assets_threadpool

def cmd_etl(args):

    ASL = Assetlist(profile=True)

    if args.ptc_coords:
        ASL.ptc_coords=True

    if args.cif_and_chains:
        ASL.cif_and_chains=True

    if args.cif_updated:
        ASL.cif_updated=True

    if args.factors_and_ligands:
        ASL.factors_and_ligands=True

    if args.png_thumbnail:
        ASL.png_thumbnail=True

    #All structures
    if args.obtain_all_structures:

        obtain_assets_threadpool(
            [],
            ASL,
            workers=16,
            get_all=True,
            overwrite=args.overwrte or False
        )

    #Single structure
    if args.rcsb_id:
        RCSB_ID = str(args.rcsb_id)
        loop    = asyncio.get_event_loop()
        loop.run_until_complete(
            obtain_assets(
                RCSB_ID,
                ASL,
                args.overwrite or False
            )
        )





