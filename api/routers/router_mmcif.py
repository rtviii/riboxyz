import os
import re
from django.http import HttpResponse
from ninja import Router
from ribctl.asset_manager.assets_structure import StructureAssets
from ribctl.ribosome_ops import RibosomeOps

router_mmcif = Router()
tag = "Crystallographic Data (mmCIF)"


@router_mmcif.get("/polymer", tags=[tag])
def polymer(request, rcsb_id: str, auth_asym_id: str):
    rcsb_id = rcsb_id.upper()
    RO = RibosomeOps(rcsb_id)
    chain_fullpath = os.path.join(
        RO.assets.paths.chains_dir, f"{rcsb_id}_{auth_asym_id}.cif"
    )

    if not os.path.exists(chain_fullpath):
        return HttpResponse(status=404)

    with open(chain_fullpath, "rb") as file:
        file_content = file.read()

    response = HttpResponse(file_content, content_type="chemical/x-mmcif")
    response["Content-Disposition"] = (
        f'attachment; filename="{os.path.basename(chain_fullpath)}"'
    )

    return response


@router_mmcif.get("/nonpolymer", tags=[tag])
def nonpolymer(request, rcsb_id: str, chemicalId: str):
    rcsb_id = rcsb_id.upper()
    p = StructureAssets(rcsb_id).paths.nonpoly_entity(chemicalId.upper())

    if not os.path.exists(p):
        return HttpResponse(status=404)

    with open(p, "rb") as file:
        file_content = file.read()

    response = HttpResponse(file_content, content_type="chemical/x-mmcif")
    response["Content-Disposition"] = f'attachment; filename="{os.path.basename(p)}"'

    return response
