# ! This should possibly be renamed as something like "ops" or "utils". Keep ligand specific for now 
import os
from pprint import pprint
from typing import Optional
import click
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.mod_extract_bsites import bsite_ligand, BindingSite, BindingSiteChain
ce = click.echo

@click.command()
@click.pass_context
@click.argument("chem_id", required=True, type=str)
@click.argument("radius",  type=float, required=False)
def lig(ctx, chem_id, radius):
    rcsb_id = ctx.obj['rcsb_id']
    print(bsite_ligand(chem_id, rcsb_id, radius).model_dump_json())