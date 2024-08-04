# ! This should possibly be renamed as something like "ops" or "utils". Keep ligand specific for now 
import json
import os
from pprint import pprint
import click
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.libbsite import bsite_ligand, BindingSite, BindingSiteChain, bsite_transpose
ce = click.echo


@click.group()
@click.pass_context
def lig(ctx):
    pass

@lig.command()
@click.argument("rcsb_id", required=True, type=str)
@click.argument("chem_id", required=True, type=str)
@click.argument("radius",  type=float, required=False, default=5)
@click.option("--save", is_flag=True,default=False)
def nbhd(rcsb_id, chem_id, radius, save):
    bsite = bsite_ligand(chem_id, rcsb_id, radius)
    if save:
        with open(RibosomeOps(rcsb_id).paths.binding_site(chem_id), 'w') as outfile:
            json.dump(bsite.model_dump(), outfile, indent=4)
            ce("Saved: {}".format( RibosomeOps(rcsb_id).paths.binding_site(chem_id)))

        

    
@lig.command()
@click.pass_context
@click.argument("lig_chem_id", required=True , type=str)
@click.argument("source_struct", required=True , type=str)
@click.argument("target_struct", required=True , type=str)
def transpose(ctx, lig_chem_id, source_struct, target_struct):
    ligpath = RibosomeOps(source_struct.upper()).paths.binding_site(lig_chem_id.upper())
    with open(ligpath, 'r') as ligfile:
        data = json.load(ligfile)
        lig = data

    target     = RibosomeOps(target_struct).profile()
    transposed = bsite_transpose(source_struct,target_struct, lig).model_dump()

    pprint(transposed)