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
@click.argument("chem_id", required=True, type=str)
@click.argument("rcsb_id", required=True, type=str)
@click.argument("radius",  type=float, required=False, default=10)
@click.option("--save", is_flag=True,default=False)
def nbhd(chem_id, rcsb_id,  radius, save):
    bsite = bsite_ligand(chem_id, rcsb_id, radius)
    if save:
        with open(RibosomeOps(rcsb_id).paths.binding_site(chem_id), 'w') as outfile:
            json.dump(bsite.model_dump(), outfile, indent=4)
            ce("Saved: {}".format(RibosomeOps(rcsb_id).paths.binding_site(chem_id)))

    pprint(bsite.model_dump())

        

    
@lig.command()
@click.pass_context
@click.argument("chem_id", required=True , type=str)
@click.argument("source_struct", required=True , type=str)
@click.argument("target_struct", required=True , type=str)
@click.argument("radius", required=True , type=float)
@click.option("--save", is_flag=True,default=False)
def transpose(ctx, chem_id, source_struct, target_struct, radius, save):

    bsite      = bsite_ligand(chem_id, source_struct, radius)
    transposed = bsite_transpose(source_struct,target_struct, bsite).model_dump()

    if save:
        transposed_path = RibosomeOps(target_struct).paths.binding_site_prediction(chem_id, source_struct)
        with open(transposed_path, 'w') as outfile:
            json.dump(transposed, outfile, indent=4)
            ce("Saved: {}".format(transposed_path))
    return transposed
