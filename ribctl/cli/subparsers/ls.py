import os
import click
from ribctl import RIBETL_DATA
from ribctl.etl.assets_structure import StructureAssets, RibosomeOps, Structure
from ribctl.lib.libtax import Taxid
ce = click.echo

# @click.command()
# @click.pass_context
# def ls(ctx):
#     rcsb_id = ctx.obj['rcsb_id']
#     if rcsb_id != None:
#         if "." in rcsb_id:
#             rcsb_id, auth_asym_id = rcsb_id.split(".")
#             RA    = RibosomeOps(rcsb_id)
#             chain = RA.get_poly_by_auth_asym_id(auth_asym_id)
#             if chain != None:
#                 ce(chain.model_dump_json())
#             else:
#                 ce(f"Chain {auth_asym_id} not found in {rcsb_id}")
#         else:
#             ce(RibosomeOps(rcsb_id).profile().model_dump_json())
#     else:
#         all_structs = os.listdir(RIBETL_DATA)
#         ce(all_structs)

@click.command()
@click.option('-t', '--taxid', type=int, help='List structures for descendants of this taxid')
@click.pass_context
def ls(ctx, taxid):
    rcsb_id = None

    if 'rcsb_id' in ctx.obj:
        rcsb_id = ctx.obj['rcsb_id'].upper()
    
    if rcsb_id is not None:
        if "." in rcsb_id:
            rcsb_id, auth_asym_id = rcsb_id.split(".")
            RA = RibosomeOps(rcsb_id)
            chain = RA.get_poly_by_auth_asym_id(auth_asym_id)
            if chain is not None:
                ce(chain.model_dump_json())
            else:
                ce(f"Chain {auth_asym_id} not found in {rcsb_id}")
        else:
            ce(RibosomeOps(rcsb_id).profile().model_dump_json())

    elif taxid is not None:
        all_structs = os.listdir(RIBETL_DATA)
        pdbid_taxid_tuples = []

        for struct in all_structs:
            try:
                ribosome_Assets = RibosomeOps(struct).profile()
                pdbid_taxid_tuples.append((ribosome_Assets.rcsb_id, ribosome_Assets.src_organism_ids[0]))
            except:
                continue

        descendants = [
            struct for struct, struct_taxid in pdbid_taxid_tuples
            if Taxid.is_descendant_of(taxid, struct_taxid)
        ]

        for struct in sorted(descendants):
            ce(struct)
    else:
        all_structs = os.listdir(RIBETL_DATA)
        for struct in sorted(all_structs):
            ce(struct)

# ... existing code ...
