import datetime
import json
import os
from pprint import pprint
import typing
from django.http import  JsonResponse, HttpResponseServerError
from ninja import Router, Schema
from pydantic import BaseModel
from neo4j_ribosome.db_lib_reader import dbqueries
from ribctl import ASSETS, ASSETS_PATH
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.info import StructureCompositionStats, run_composition_stats
from ribctl.lib.schema.types_ribosome import  CytosolicProteinClass, CytosolicRNAClass, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, MitochondrialProteinClass, MitochondrialRNAClass, PTCInfo, Polymer, PolymerClass, PolymerClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure, tRNA
from ribctl.lib.libtax import Taxid 

structure_router = Router()
TAG              = "structures"

@structure_router.get("/all_rcsb_ids", response=list[str], tags=[TAG])
def all_rcsb_ids(request):
    return dbqueries.all_ids()

@structure_router.get("/polymer_classes_stats", response=list[tuple[str,int]], tags=[TAG])
def polymer_classes_stats(request):
    return dbqueries.polymer_classes_stats()

@structure_router.get("/tax_dict", response=dict, tags=[TAG])
def tax_dict(request):
    """Returns a dictionary of all taxonomic IDs present in the database as keys, [ scientific name, corresponding superkingdom ] as the values."""
    td = dbqueries.tax_dict()
    _  = {}
    for kvp in td:
        _[kvp[0]] = kvp[1]
    return _

@structure_router.get("/polymer_classification_report",  tags=[TAG])
def polymer_classification_report(request, rcsb_id:str, auth_asym_id:str):
    """Returns a dictionary of all taxonomic IDs present in the database as keys, [ scientific name, corresponding superkingdom ] as the values."""
    if os.path.exists(RibosomeOps(rcsb_id).paths.classification_report):
        with open(RibosomeOps(rcsb_id).paths.classification_report, 'r') as f:
            return json.load(f)[auth_asym_id]
    else :
        return []

@structure_router.get("/structure_composition_stats", response=StructureCompositionStats, tags=[TAG])
def structure_composition_stats(request):
    filename = os.path.join(ASSETS_PATH, "structure_composition_stats.json")
    if os.path.exists(filename):
        creation_time = os.path.getctime(filename)
        creation_date = datetime.datetime.fromtimestamp(creation_time)
        #TODO:  if past the last update date, rerun
    else:
        run_composition_stats()

    with open(filename, 'r') as infile:
        return json.load(infile)

@structure_router.get("/random_profile", response=RibosomeStructure, tags=[TAG])
def random_profile(request):
    return RibosomeStructure.model_validate(dbqueries.random_structure()[0])

@structure_router.get('/list_polymers_filtered_by_polymer_class', response=dict,  tags=[TAG])
def polymers_by_polymer_class(request,
      polymer_class: PolymerClass,
      page  = 1,
      ):

    polymers, count = dbqueries.list_polymers_filtered_by_polymer_class(page, polymer_class)[0]
    return { "polymers":polymers, "count": count }

@structure_router.get('/list_polymers_by_structure', response=dict,  tags=[TAG])
def polymers_by_structure(request,
      page  = 1,
      search          = None,
      year            = None,
      resolution      = None,
      polymer_classes = None,
      source_taxa     = None,
      host_taxa       = None):

                
    def parse_empty_or_int(_:str):
        if _ != '':
            return int(_)
        else:
            return None

    def parse_empty_or_float(_:str):
        if _ != '':
            return float(_)
        else:
            return None


    year            = None if year            == "" else list(map(parse_empty_or_int,year.split(",")))
    resolution      = None if resolution      == "" else list(map(parse_empty_or_float,resolution.split(",")))
    host_taxa       = None if host_taxa       == "" else list(map(parse_empty_or_int,host_taxa.split(","))) if host_taxa else None
    source_taxa     = None if source_taxa     == "" else list(map(parse_empty_or_int,source_taxa.split(","))) if source_taxa else None
    polymer_classes = None if polymer_classes == "" else list(map(lambda _: PolymerClass(_), polymer_classes.split(",")))

    qreturn =  dbqueries.list_polymers_filtered_by_structure(page, search, year, resolution, polymer_classes, source_taxa, host_taxa)

    if len(qreturn) < 1:
        print("Found none. returning empty", { "polymers":[], "count": 0 })
        return { "polymers":[], "count": 0 }
    else:
        polymers, count = qreturn[0]
        print("RETURNING POLYMERS actual len", len( polymers ), count)
        return { "polymers":polymers, "count": count }

@structure_router.get('/ptc',  tags=[TAG], response=PTCInfo)
def ptc(request, rcsb_id:str):
    rcsb_id = str.upper(rcsb_id)
    return RibosomeOps(rcsb_id).ptc()

@structure_router.get('/list_ligands',response=list[tuple[dict,list[dict]]] , tags=[TAG])
def list_ligands(request):
    return dbqueries.list_ligands()

@structure_router.get('/list', response=dict,  tags=[TAG])
def filter_list(request,
      page            = 1,
      search          = None,
      year            = None,
      resolution      = None,
      polymer_classes = None,
      source_taxa     = None,
      host_taxa       = None):

    def parse_empty_or_int(_:str):
        if _ != '':
            return int(_)
        else:
            return None

    def parse_empty_or_float(_:str):
        if _ != '':
            return float(_)
        else:
            return None


    year            = None if year            == "" else list(map(parse_empty_or_int,year.split(",")))
    resolution      = None if resolution      == "" else list(map(parse_empty_or_float,resolution.split(",")))
    host_taxa       = None if host_taxa       == "" else list(map(parse_empty_or_int,host_taxa.split(","))) if host_taxa else None
    source_taxa     = None if source_taxa     == "" else list(map(parse_empty_or_int,source_taxa.split(","))) if source_taxa else None
    polymer_classes = None if polymer_classes == "" else list(map(lambda _: PolymerClass(_), polymer_classes.split(",")))

    structures, count = dbqueries.list_structs_filtered(int(page), search, year, resolution, polymer_classes, source_taxa, host_taxa)[0]
    structures_validated = list(map(lambda r: RibosomeStructure.model_validate(r), structures))
    return { "structures":structures_validated, "count": count }


@structure_router.get('/profile', response=RibosomeStructure, tags=[TAG],)
def structure_profile(request,rcsb_id:str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])

    try:
        with open(RibosomeOps(rcsb_id).paths.profile, 'r') as f:
            return JsonResponse(json.load(f))
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))

@structure_router.get('/ptc', response=dict, tags=[TAG],)
def structure_ptc(request,rcsb_id:str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    try:
        ptc = RibosomeOps(rcsb_id).ptc()
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))
     
class ChainsByStruct(Schema):
    class PolymerByStruct(Schema):
        nomenclature: list[PolymerClass]
        auth_asym_id: str
        entity_poly_polymer_type: str
        entity_poly_seq_length: int

    polymers: list[PolymerByStruct]
    rcsb_id : str
   
@structure_router.get('/chains_by_struct', response=list[ChainsByStruct], tags=[TAG])
def chains_by_struct(request):
    structs_response = dbqueries.list_chains_by_struct()
    return structs_response

class NomenclatureSet(Schema)  : 
      ElongationFactorClass    : list[str]
      InitiationFactorClass    : list[str]
      CytosolicProteinClass    : list[str]
      MitochondrialProteinClass: list[str]
      CytosolicRNAClass        : list[str]
      MitochondrialRNAClass    : list[str]
      tRNAClass                : list[str]


@structure_router.get('/list_nomenclature', response=NomenclatureSet)
def polymer_classes_nomenclature(request):
    classes = {
        "ElongationFactorClass"    : [e.value for e in ElongationFactorClass],
        "InitiationFactorClass"    : [e.value for e in InitiationFactorClass],
        "CytosolicProteinClass"    : [e.value for e in CytosolicProteinClass ],
        "MitochondrialProteinClass": [e.value for e in MitochondrialProteinClass ],
        "CytosolicRNAClass"        : [e.value for e in CytosolicRNAClass ],
        "MitochondrialRNAClass"    : [e.value for e in MitochondrialRNAClass ],
        "tRNAClass"                : [e.value for e in tRNA ],
    }
    return  classes


@structure_router.get('/list_source_taxa', response=list[dict], tags=[TAG])
def list_source_taxa(request, source_or_host:typing.Literal["source", "host"]):
    """This endpoint informs the frontend about which tax ids are present in the database.
    Used for filters/search. """

    tax_ids = dbqueries.get_taxa(source_or_host)
    def normalize_tax_list_to_dict(tax_list:list[int]):
        # TODO: This can be done a lot better with recursive descent. 
        # TODO: Error hnadling?
        # You could also parametrize the "include_only" given that all Taxid handles that.

        def inject_species(node, S:int, F:int, L:int):
            """This acts on the family node"""
            global nodes
            if node['taxid'] == F:
                if len(list(filter(lambda subnode: subnode['taxid'] == S, node['children'])) ) < 1:
                    node['children'].append({'taxid': S, "value":S, 'title': Taxid.get_name(S),  })
            return node

        def inject_families(node, S:int, F:int, K:int):
            """This acts on the superkingdom node"""
            global nodes
            if node['taxid'] == K:
                if len(list(filter(lambda subnode: subnode['taxid'] == F, node['children'])) ) < 1:
                    node['children'].append({'taxid': F, "value":F,'title': Taxid.get_name(F), 'children': []})
                list(map(lambda node: inject_species(node,S, F, K), node['children']))
            return node

        normalized_taxa = []
        for tax in tax_list:
            p = Taxid.get_lineage(tax, include_only=['superkingdom', 'family', 'species'])
            K,F,S = p
            if len( list(filter(lambda obj: obj['taxid'] == K, normalized_taxa)) ) < 1:
                normalized_taxa.append({'taxid': K, "value":K,'title': Taxid.get_name(K), "children": []})
            
            list(map(lambda node: inject_families(node,S, F, K), normalized_taxa))

        return normalized_taxa
       
    
    return normalize_tax_list_to_dict(tax_ids)
