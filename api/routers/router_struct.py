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

router_structures = Router()
TAG_STRUCTURES = "structures"
TAG_LANDMARKS = "landmarks"


@router_structures.get("/all_rcsb_ids", response=list[str], tags=[TAG_STRUCTURES])
def all_rcsb_ids(request):
    return dbqueries.all_ids()

@router_structures.get( "/polymer_classes_stats", response=list[tuple[str, int]], tags=[TAG_STRUCTURES] )
def polymer_classes_stats(request):
    return dbqueries.polymer_classes_stats()


@router_structures.get("/tax_dict", response=dict, tags=[TAG_STRUCTURES])
def tax_dict(request):
    """Returns a dictionary of all taxonomic IDs present in the database as keys, [ scientific name, corresponding superkingdom ] as the values."""
    td = dbqueries.tax_dict()
    _ = {}
    for kvp in td:
        _[kvp[0]] = kvp[1]
    return _


@router_structures.get("/polymer_classification_report", tags=[TAG_STRUCTURES])
def polymer_classification_report(request, rcsb_id: str):
    if os.path.exists(RibosomeOps(rcsb_id).assets.paths.classification_report):
        with open(RibosomeOps(rcsb_id).assets.paths.classification_report, "r") as f:
            return json.load(f)
    else:
        return []


@router_structures.get("/random_profile", response=RibosomeStructure, tags=[TAG_STRUCTURES])
def random_profile(request):
    return RibosomeStructure.model_validate(dbqueries.random_structure()[0])




@router_structures.post("/list_structures", response=dict, tags=[TAG_STRUCTURES])
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




@router_structures.get("/structures_overview", response=list[dict], tags=[TAG_STRUCTURES])
def overview(request):
    return dbqueries.structures_overview()


@router_structures.get(
    "/profile",
    response=RibosomeStructure,
    tags=[TAG_STRUCTURES],
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


@router_structures.get(
    "/ptc",
    response=PTCInfo,
    tags=[TAG_STRUCTURES, TAG_LANDMARKS ],
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


@router_structures.get(
    "/constriction_site",
    response=ConstrictionSite,
    tags=[TAG_STRUCTURES],
)
def constriction_site(request, rcsb_id: str):
    params = dict(request.GET)
    rcsb_id = str.upper(params["rcsb_id"][0])
    try:
        # constriction = AssetPathManager().load_model(
        #     rcsb_id, AssetType.CONSTRICTION_SITE
        # )
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


class ChainsByStruct(Schema):
    class PolymerByStruct(Schema):
        nomenclature: list[PolymerClass]
        auth_asym_id: str
        entity_poly_polymer_type: str
        entity_poly_seq_length: int

    polymers: list[PolymerByStruct]
    rcsb_id: str


@router_structures.get("/chains_by_struct", response=list[ChainsByStruct], tags=[TAG_STRUCTURES])
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


@router_structures.get("/list_nomenclature", response=NomenclatureSet)
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


# @router_structures.get("/list_source_taxa", response=list[dict], tags=[TAG_STRUCTURES])
# def list_source_taxa(request, source_or_host: typing.Literal["source", "host"]):
#     """This endpoint informs the frontend about which tax ids are present in the database.
#     Used for filters/search."""

#     tax_ids = dbqueries.get_taxa(source_or_host)
#     def normalize_tax_list_to_dict(tax_list: list[int]):
#         # TODO: This can be done a lot better with recursive descent.
#         # TODO: Error hnadling?
#         # You could also parametrize the "include_only" given that all Taxid handles that.

#         def inject_species(node, S: int, F: int, L: int):
#             """This acts on the family node"""
#             global nodes
#             if node["taxid"] == F:
#                 if (
#                     len(
#                         list(
#                             filter(
#                                 lambda subnode: subnode["taxid"] == S, node["children"]
#                             )
#                         )
#                     )
#                     < 1
#                 ):
#                     node["children"].append(
#                         {
#                             "taxid": S,
#                             "value": S,
#                             "title": Taxid.get_name(S),
#                         }
#                     )
#             return node

#         def inject_families(node, S: int, F: int, K: int):
#             """This acts on the superkingdom node"""
#             global nodes
#             if node["taxid"] == K:
#                 if (
#                     len(
#                         list(
#                             filter(
#                                 lambda subnode: subnode["taxid"] == F, node["children"]
#                             )
#                         )
#                     )
#                     < 1
#                 ):
#                     node["children"].append(
#                         {
#                             "taxid": F,
#                             "value": F,
#                             "title": Taxid.get_name(F),
#                             "children": [],
#                         }
#                     )
#                 list(map(lambda node: inject_species(node, S, F, K), node["children"]))
#             return node

#         normalized_taxa = []
#         for tax in tax_list:
#             p = Taxid.get_lineage(
#                 tax, include_only=["superkingdom", "family", "species"]
#             )
#             K, F, S = p
#             if len(list(filter(lambda obj: obj["taxid"] == K, normalized_taxa))) < 1:
#                 normalized_taxa.append(
#                     {"taxid": K, "value": K, "title": Taxid.get_name(K), "children": []}
#                 )
# # 
#             list(map(lambda node: inject_families(node, S, F, K), normalized_taxa))

#         return normalized_taxa

#     return normalize_tax_list_to_dict(tax_ids)

@router_structures.get("/list_source_taxa", response=list[dict], tags=[TAG_STRUCTURES])
def list_source_taxa(request, source_or_host: typing.Literal["source", "host"]):
    """This endpoint informs the frontend about which tax ids are present in the database.
    Used for filters/search."""

    def create_node(taxid: int, include_children: bool = True) -> dict:
        """Create a tree node with standard format"""
        node = {
            "taxid": taxid,
            "value": taxid,
            "title": Taxid.get_name(taxid)
        }
        if include_children:
            node["children"] = []
        return node

    def build_taxonomy_tree(tax_ids: list[int]) -> list[dict]:
        # Dictionary to store our hierarchy
        taxonomy_tree = {}
        
        for taxid in tax_ids:
            try:
                # Get lineage for current taxid
                lineage = Taxid.get_lineage(
                    taxid, 
                    include_only=["superkingdom", "family", "species"]
                )
                
                if len(lineage) != 3:
                    continue  # Skip incomplete lineages
                    
                superkingdom_id, family_id, species_id = lineage
                
                # Add superkingdom if not exists
                if superkingdom_id not in taxonomy_tree:
                    taxonomy_tree[superkingdom_id] = create_node(superkingdom_id)
                
                # Add family if not exists
                superkingdom_node = taxonomy_tree[superkingdom_id]
                family_exists = any(child["taxid"] == family_id 
                                  for child in superkingdom_node["children"])
                
                if not family_exists:
                    superkingdom_node["children"].append(create_node(family_id))
                
                # Add species if not exists
                family_node = next(child for child in superkingdom_node["children"] 
                                 if child["taxid"] == family_id)
                
                species_exists = any(child["taxid"] == species_id 
                                   for child in family_node.get("children", []))
                
                if not species_exists:
                    if "children" not in family_node:
                        family_node["children"] = []
                    family_node["children"].append(create_node(species_id, False))
                    
            except Exception as e:
                print(f"Error processing taxid {taxid}: {str(e)}")
                continue
        
        return list(taxonomy_tree.values())

    tax_ids = dbqueries.get_taxa(source_or_host)
    return build_taxonomy_tree(tax_ids)


@router_structures.get("/tunnel_geometry")
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


@router_structures.get("/cylinder_residues")
def cylinder_residues(request):
    try:
        with open("/home/rtviii/dev/npet-cg-sim/cylinder_residues.json", "rb") as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_structures.get("/half_cylinder_residues")
def half_cylinder_residues(request):
    try:
        with open(
            "/home/rtviii/dev/npet-cg-sim/half_cylinder_residues.json", "rb"
        ) as f:
            map = json.load(f)
        return Response(map)
    except IOError:
        print("Couldn not read file")


@router_structures.get("/landmarks/helices/{rcsb_id}")
def get_helices(request, rcsb_id: str):
    rcsb_id = rcsb_id.upper()
    file_path = os.path.join("/home/rtviii/dev/riboxyz/7K00_rrna_helices.json")
    try:
        with open(file_path, "rb") as f:
            result = json.load(f)
        return Response(result)
    except IOError as e:
        return Response({"error": f"Error reading the helices file {e}"}, status=500)
