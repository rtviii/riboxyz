import os
import click
from ribctl import RIBETL_DATA
from ribctl.etl.etl_ribosome_ops import RibosomeOps, Structure
ce = click.echo

@click.command()
@click.pass_context
def ls(ctx):

    rcsb_id = ctx.obj['rcsb_id']
    if rcsb_id != None:
        if "." in rcsb_id:
            rcsb_id, auth_asym_id = rcsb_id.split(".")
            RA    = RibosomeOps(rcsb_id)
            chain = RA.get_poly_by_auth_asym_id(auth_asym_id)
            if chain != None:
                ce(chain.model_dump_json())
            else:
                ce(f"Chain {auth_asym_id} not found in {rcsb_id}")
        else:
            ce(RibosomeOps(rcsb_id).profile().model_dump_json())
    else:
        all_structs = os.listdir(RIBETL_DATA)
        ce(all_structs)


    # print(all_structs)




# if args.struct != None: