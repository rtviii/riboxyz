import datetime
import json
import pickle
from urllib.request import Request
from ninja import Router, Body
import os
from pprint import pprint
import typing
from django.http import JsonResponse, HttpResponseServerError
from ninja import Path, Router, Schema
import pandas
from pydantic import BaseModel, Field, ValidationError
from typing import Optional, List
from pydantic import BaseModel
from neo4j_ribosome.db_lib_reader import (
    PolymersFilterParams,
    StructureFilterParams,
    dbqueries,
)
from ribctl import ASSETS, ASSETS_PATH, RIBETL_DATA
from ribctl.asset_manager.asset_manager import AssetPathManager
from ribctl.asset_manager.types import AssetType
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.info import StructureCompositionStats, run_composition_stats
from ribctl.lib.types.polymer import (
    CytosolicProteinClass,
    CytosolicRNAClass,
    ElongationFactorClass,
    InitiationFactorClass,
    LifecycleFactorClass,
    MitochondrialProteinClass,
    MitochondrialRNAClass,
    PolymerClass,
    PolynucleotideClass,
    PolynucleotideClass,
    PolypeptideClass,
    ProteinClass,
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

structure_router = Router()
TAG = "structures"


@structure_router.get("/all_rcsb_ids", response=list[str], tags=[TAG])
def all_rcsb_ids(request):
    return dbqueries.all_ids()


@structure_router.get(
    "/polymer_classes_stats", response=list[tuple[str, int]], tags=[TAG]
)
def polymer_classes_stats(request):
    return dbqueries.polymer_classes_stats()


@structure_router.get("/tax_dict", response=dict, tags=[TAG])
def tax_dict(request):
    """Returns a dictionary of all taxonomic IDs present in the database as keys, [ scientific name, corresponding superkingdom ] as the values."""
    td = dbqueries.tax_dict()
    _ = {}
    for kvp in td:
        _[kvp[0]] = kvp[1]
    return _


@structure_router.get("/polymer_classification_report", tags=[TAG])
def polymer_classification_report(request, rcsb_id: str):
    if os.path.exists(RibosomeOps(rcsb_id).assets.paths.classification_report):
        with open(RibosomeOps(rcsb_id).assets.paths.classification_report, "r") as f:
            return json.load(f)
    else:
        return []


@structure_router.get(
    "/structure_composition_stats", response=StructureCompositionStats, tags=[TAG]
)
def structure_composition_stats(request):
    filename = os.path.join(ASSETS_PATH, "structure_composition_stats.json")
    if os.path.exists(filename):
        creation_time = os.path.getctime(filename)
        creation_date = datetime.datetime.fromtimestamp(creation_time)
        # TODO:  if past the last update date, rerun
    else:
        run_composition_stats()

    with open(filename, "r") as infile:
        return json.load(infile)


@structure_router.get("/random_profile", response=RibosomeStructure, tags=[TAG])
def random_profile(request):
    return RibosomeStructure.model_validate(dbqueries.random_structure()[0])


@structure_router.get(
    "/list_ligands", response=list[tuple[dict, list[dict]]], tags=[TAG]
)
def list_ligands(request):
    return dbqueries.list_ligands()


@structure_router.post("/list_structures", response=dict, tags=[TAG])
def list_structures(request, filters: StructureFilterParams):
    parsed_filters = StructureFilterParams(**json.loads(request.body))
    print(parsed_filters)
    structures, next_cursor, total_count = dbqueries.list_structs_filtered(
        parsed_filters
    )
    structures_validated = [
        RibosomeStructureMetadata.model_validate(s) for s in structures
    ]
    return {
        "structures": structures_validated,
        "next_cursor": next_cursor,
        "total_count": total_count,
    }


@structure_router.post("/list_polymers", response=dict, tags=[TAG])
def list_polymers(request, filters: PolymersFilterParams):
    parsed_filters = PolymersFilterParams(**json.loads(request.body))
    print(parsed_filters)
    polymers, next_cursor, total_polymers_count, total_structures_count = (
        dbqueries.list_polymers_filtered(parsed_filters)
    )
    polymers_validated = [Polymer.model_validate(p) for p in polymers]
    return {
        "polymers": polymers_validated,
        "next_cursor": next_cursor,
        "total_polymers_count": total_polymers_count,
        "total_structures_count": total_structures_count,
    }


@structure_router.get("/structures_overview", response=list[dict], tags=[TAG])
def overview(request):
    return dbqueries.structures_overview()


@structure_router.get(
    "/profile",
    response=RibosomeStructure,
    tags=[TAG],
)
def structure_profile(request, rcsb_id: str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])

    try:
        with open(RibosomeOps(rcsb_id).assets.paths.profile, "r") as f:
            return JsonResponse(json.load(f))
    except Exception as e:
        return HttpResponseServerError(
            "Failed to find structure profile {}:\n\n{}".format(rcsb_id, e)
        )


@structure_router.get(
    "/ptc",
    response=PTCInfo,
    tags=[TAG],
)
def structure_ptc(request, rcsb_id: str):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        ptc = AssetPathManager().load_model(rcsb_id, AssetType.PTC)
        if not ptc:
            return JsonResponse({"error": "No PTC found for {}".format(rcsb_id)})
        return JsonResponse(ptc.model_dump())
    except Exception as e:
        return HttpResponseServerError(e)


@structure_router.get(
    "/constriction_site",
    response=ConstrictionSite,
    tags=[TAG],
)
def constriction_site(request, rcsb_id: str):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        constriction = AssetPathManager().load_model(
            rcsb_id, AssetType.CONSTRICTION_SITE
        )
        if not constriction:
            return JsonResponse(
                {"error": "No constriction site found for {}".format(rcsb_id)}
            )
        return JsonResponse(constriction.model_dump())
    except Exception as e:
        return HttpResponseServerError(e)


class ChainsByStruct(Schema):
    class PolymerByStruct(Schema):
        nomenclature: list[PolymerClass]
        auth_asym_id: str
        entity_poly_polymer_type: str
        entity_poly_seq_length: int

    polymers: list[PolymerByStruct]
    rcsb_id: str


@structure_router.get("/chains_by_struct", response=list[ChainsByStruct], tags=[TAG])
def chains_by_struct(request):
    structs_response = dbqueries.list_chains_by_struct()
    return structs_response


class NomenclatureSet(Schema):
    ElongationFactorClass: list[str]
    InitiationFactorClass: list[str]
    CytosolicProteinClass: list[str]
    MitochondrialProteinClass: list[str]
    CytosolicRNAClass: list[str]
    MitochondrialRNAClass: list[str]
    tRNAClass: list[str]


@structure_router.get("/list_nomenclature", response=NomenclatureSet)
def polymer_classes_nomenclature(request):
    classes = {
        "ElongationFactorClass": [e.value for e in ElongationFactorClass],
        "InitiationFactorClass": [e.value for e in InitiationFactorClass],
        "CytosolicProteinClass": [e.value for e in CytosolicProteinClass],
        "MitochondrialProteinClass": [e.value for e in MitochondrialProteinClass],
        "CytosolicRNAClass": [e.value for e in CytosolicRNAClass],
        "MitochondrialRNAClass": [e.value for e in MitochondrialRNAClass],
        "tRNAClass": [e.value for e in tRNA],
    }
    return classes


@structure_router.get("/list_source_taxa", response=list[dict], tags=[TAG])
def list_source_taxa(request, source_or_host: typing.Literal["source", "host"]):
    """This endpoint informs the frontend about which tax ids are present in the database.
    Used for filters/search."""

    tax_ids = dbqueries.get_taxa(source_or_host)

    def normalize_tax_list_to_dict(tax_list: list[int]):
        # TODO: This can be done a lot better with recursive descent.
        # TODO: Error hnadling?
        # You could also parametrize the "include_only" given that all Taxid handles that.

        def inject_species(node, S: int, F: int, L: int):
            """This acts on the family node"""
            global nodes
            if node["taxid"] == F:
                if (
                    len(
                        list(
                            filter(
                                lambda subnode: subnode["taxid"] == S, node["children"]
                            )
                        )
                    )
                    < 1
                ):
                    node["children"].append(
                        {
                            "taxid": S,
                            "value": S,
                            "title": Taxid.get_name(S),
                        }
                    )
            return node

        def inject_families(node, S: int, F: int, K: int):
            """This acts on the superkingdom node"""
            global nodes
            if node["taxid"] == K:
                if (
                    len(
                        list(
                            filter(
                                lambda subnode: subnode["taxid"] == F, node["children"]
                            )
                        )
                    )
                    < 1
                ):
                    node["children"].append(
                        {
                            "taxid": F,
                            "value": F,
                            "title": Taxid.get_name(F),
                            "children": [],
                        }
                    )
                list(map(lambda node: inject_species(node, S, F, K), node["children"]))
            return node

        normalized_taxa = []
        for tax in tax_list:
            p = Taxid.get_lineage(
                tax, include_only=["superkingdom", "family", "species"]
            )
            K, F, S = p
            if len(list(filter(lambda obj: obj["taxid"] == K, normalized_taxa))) < 1:
                normalized_taxa.append(
                    {"taxid": K, "value": K, "title": Taxid.get_name(K), "children": []}
                )

            list(map(lambda node: inject_families(node, S, F, K), normalized_taxa))

        return normalized_taxa

    return normalize_tax_list_to_dict(tax_ids)


from django.http import FileResponse
from ninja.responses import Response


@structure_router.get("/tunnel_geometry")
def get_shape(request, rcsb_id: str, is_ascii: bool = False):
    rcsb_id = rcsb_id.upper()
    filename = (
        "{}_poisson_recon.ply".format(rcsb_id)
        if not is_ascii
        else "{}_poisson_recon_ascii.ply".format(rcsb_id)
    )
    file_path = os.path.join(RIBETL_DATA, rcsb_id, "TUNNELS", filename)

    if not os.path.exists(file_path):
        return Response({"error": "Shape file not found"}, status=404)
    try:
        file = open(file_path, "rb")
        return FileResponse(
            file, content_type="application/octet-stream", filename=filename
        )
    except IOError:
        return Response({"error": "Error reading the shape file"}, status=500)


@structure_router.get("/cylinder_residues")
def cylinder_residues(request):
    try:
        with open("/home/rtviii/dev/npet-cg-sim/cylinder_residues.json", "rb") as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@structure_router.get("/half_cylinder_residues")
def half_cylinder_residues(request):
    try:
        with open(
            "/home/rtviii/dev/npet-cg-sim/half_cylinder_residues.json", "rb"
        ) as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@structure_router.get("/landmarks/helices/{rcsb_id}")
def get_helices(request, rcsb_id: str):
    rcsb_id = rcsb_id.upper()
    file_path = os.path.join("/home/rtviii/dev/riboxyz/7K00_rrna_helices.json")
    try:
        with open(file_path, "rb") as f:
            result = json.load(f)
        return Response(result)
    except IOError as e:
        return Response({"error": f"Error reading the helices file {e}"}, status=500)
