import datetime
import json
from urllib.request import Request
from ninja import Router, Body
import os
from pprint import pprint
import typing
from django.http import JsonResponse, HttpResponseServerError
from ninja import Path, Router, Schema
from neo4j_ribosome.db_lib_reader import (
    PolymersFilterParams,
    StructureFilterParams,
    dbqueries,
)

from django.http import FileResponse
from ninja.responses import Response
from ribctl import ASSETS, ASSETS_PATH, RIBETL_DATA
from ribctl.asset_manager.asset_manager import RibosomeAssetManager
from ribctl.asset_manager.asset_registry import AssetRegistry
from ribctl.asset_manager.asset_types import AssetType
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.info import StructureCompositionStats, run_composition_stats
from ribctl.lib.types.polymer import (
    CytosolicProteinClass,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    MitochondrialProteinClass,
    MitochondrialRNAClass,
    PolymerClass,
    tRNA,
)
from ribctl.lib.schema.types_ribosome import (
    ConstrictionSite,
    PTCInfo,
    RibosomeStructure,
    RibosomeStructureMetadata,
    Polymer,
    RibosomeStructure,
    RibosomeStructureMetadata,
    RibosomeStructureMetadata,
)
from ribctl.lib.libtax import Taxid

router_loci = Router()
TAG_LOCI = "Biologically Relevant Loci & Landmarks"


@router_loci.get("/tunnel_geometry", tags=[TAG_LOCI])
def get_shape(request, rcsb_id: str, is_ascii: bool = False):
    rcsb_id = rcsb_id.upper()
    filename = (
        "{}_poisson_recon.ply".format(rcsb_id)
        if not is_ascii
        else "{}_poisson_recon_ascii.ply".format(rcsb_id)
    )
    file_path = AssetType.NPET_MESH_ASCII.get_path(rcsb_id)

    if not os.path.exists(file_path):
        return Response({"error": "Shape file not found"}, status=404)
    try:
        file = open(file_path, "rb")
        return FileResponse(
            file, content_type="application/octet-stream", filename=filename
        )
    except IOError:
        return Response({"error": "Error reading the shape file"}, status=500)


@router_loci.get("/cylinder_residues", tags=[TAG_LOCI], include_in_schema=False)
def cylinder_residues(request):
    try:
        with open("/home/rtviii/dev/npet-cg-sim/cylinder_residues.json", "rb") as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_loci.get("/half_cylinder_residues", tags=[TAG_LOCI], include_in_schema=False)
def half_cylinder_residues(request):
    try:
        with open(
            "/home/rtviii/dev/npet-cg-sim/half_cylinder_residues.json", "rb"
        ) as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_loci.get("/helices/{rcsb_id}", tags=[TAG_LOCI])
def get_helices(request, rcsb_id: str):
    rcsb_id = rcsb_id.upper()
    file_path = os.path.join("/home/rtviii/dev/riboxyz/7K00_rrna_helices.json")
    try:
        with open(file_path, "rb") as f:
            result = json.load(f)
        return Response(result)
    except IOError as e:
        return Response({"error": f"Error reading the helices file {e}"}, status=500)


@router_loci.get(
    "/ptc",
    response=PTCInfo,
    tags=[TAG_LOCI],
)
def structure_ptc(request, rcsb_id: str):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        ptc = RibosomeAssetManager().load_model(rcsb_id, AssetType.PTC)
        path = AssetType.PTC.get_path(rcsb_id)
        if not ptc:
            return JsonResponse({"error": "No PTC found for {}".format(rcsb_id)})
        return JsonResponse(ptc.model_dump())
    except Exception as e:
        return HttpResponseServerError(e)


@router_loci.get(
    "/constriction_site",
    response=ConstrictionSite,
    tags=[TAG_LOCI],
)
def constriction_site(request, rcsb_id: str,):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        constriction = RibosomeAssetManager().load_model(
            rcsb_id, AssetType.CONSTRICTION_SITE
        )
        if not constriction:
            return JsonResponse(
                {"error": "No constriction site found for {}".format(rcsb_id)}
            )
        return JsonResponse(constriction.model_dump())
    except Exception as e:
        return HttpResponseServerError(e)
