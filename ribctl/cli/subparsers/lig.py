# ! This should possibly be renamed as something like "ops" or "utils". Keep ligand specific for now 
import json
import os
from pprint import pprint
from typing import Optional
import click
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.mod_extract_bsites import bsite_ligand, BindingSite, BindingSiteChain, init_transpose_ligand
ce = click.echo



@click.group()
@click.pass_context
def lig(ctx):
    pass


@lig.command()
@click.argument("rcsb_id", required=True, type=str)
@click.argument("chem_id", required=True, type=str)
@click.argument("radius",  type=float, required=False, default=5)
def nbhd(rcsb_id, chem_id, radius):
    print(bsite_ligand(chem_id, rcsb_id, radius).model_dump_json())
    

@lig.command()
@click.pass_context
@click.argument("lig_chem_id", required=True , type=str)
@click.argument("source_struct", required=True , type=str)
@click.argument("target_struct", required=True , type=str)
def transpose(ctx, lig_chem_id, source_struct, target_struct):
    ligpath = RibosomeOps(source_struct.upper()).paths.binding_site(lig_chem_id.upper())
    with open(ligpath, 'r') as ligfile:
        data = json.load(ligfile)
        lig = BindingSite.model_validate(data)

    target     = RibosomeOps(target_struct).profile()
    transposed = init_transpose_ligand(target, lig).model_dump()

    pprint(transposed)