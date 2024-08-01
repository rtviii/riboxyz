# ! This should possibly be renamed as something like "ops" or "utils". Keep ligand specific for now 
import os
from pprint import pprint
from typing import Optional
import click
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.mod_extract_bsites import bsite_ligand, BindingSite, BindingSiteChain
from ribctl.lib.mod_transpose_bsites import init_transpose_ligand
ce = click.echo



@click.group()
@click.pass_context
def lig(ctx):
    pass


@lig.command()
@click.argument("chem_id", required=True, type=str)
@click.argument("radius",  type=float, required=False)
def nbhd(ctx, chem_id, radius):
    rcsb_id = ctx.obj['rcsb_id']
    print(bsite_ligand(chem_id, rcsb_id, radius).model_dump_json())
    


@lig.command()
@click.pass_context
@click.argument("lig_chem_id", required=True , type=str)
@click.argument("source_struct", required=True , type=str)
@click.argument("target_struct", required=True , type=str)
def transpose(ctx, target):
    

    init_transpose_ligand()