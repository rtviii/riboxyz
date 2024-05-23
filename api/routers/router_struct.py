import json
from pprint import pprint
import typing
from typing import Optional, Tuple
from django.http import  JsonResponse, HttpResponseServerError
from ninja import Router, Schema
from neo4j_ribosome.db_reader import dbqueries
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import  CytosolicProteinClass, CytosolicRNAClass, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, MitochondrialProteinClass, MitochondrialRNAClass, PolymerClass, PolynucleotideClass, PolypeptideClass, ProteinClass, RibosomeStructure, tRNA
from ribctl.lib.libtax import Taxid, ncbi

structure_router = Router()
TAG              = "structures"

@structure_router.get('/profile', response=RibosomeStructure, tags=[TAG],)
def structure_profile(request,rcsb_id:str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])

    try:
        with open(RibosomeAssets(rcsb_id)._json_profile_filepath(), 'r') as f:
            return JsonResponse(json.load(f))
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))

@structure_router.get('/ptc', response=dict, tags=[TAG],)
def structure_ptc(request,rcsb_id:str):
    """Return a `.json` profile of the given RCSB_ID structure."""

    params      = dict(request.GET)
    rcsb_id     = str.upper(params['rcsb_id'][0])
    try:
        ptc = RibosomeAssets(rcsb_id)._ptc_residues()
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))



# class StructsListSchema(Schema): 
#       search                   : Optional[str                                           ]= None
#       year                     : Optional[typing.Tuple[int | None , int | None]         ]= None
#       resolution               : Optional[typing.Tuple[float | None , float | None]     ]= None
#       polymer_classes          : Optional[list[PolynucleotideClass | PolypeptideClass ] ]= None
#       source_taxa              : Optional[list[int]                                     ]= None
#       host_taxa                : Optional[list[int]                                     ]= None

@structure_router.get('/list', response=dict,  tags=[TAG])
def filter_list(request,
      search          = None,
      year            = None,
      resolution      = None,
      polymer_classes = None,
      source_taxa     = None,
      host_taxa       = None):

    print("Got params:" )
    print("Search:", search)
    print("year:", year)
    print("resolution:", resolution)
    print("polymer_classes:", polymer_classes)
    print("source_taxa:", source_taxa)
    print("host_taxa:", host_taxa)
                
    def parse_empty_or_int(_:str):
        if _ != '':
            return int(_)
        else:
            return None



    year            = list(map(parse_empty_or_int,year.split(","))) if year else None
    resolution      = list(map(parse_empty_or_int,resolution.split(",")))  if resolution else None
    host_taxa       = list(map(parse_empty_or_int,host_taxa.split(","))) if host_taxa else None
    source_taxa     = list(map(parse_empty_or_int,source_taxa.split(","))) if source_taxa else None
    polymer_classes = list(map(PolymerClass,polymer_classes.split(","))) if polymer_classes else None

    structures, count = dbqueries.list_structs_filtered(search, year, resolution, polymer_classes, source_taxa, host_taxa)[0]
    structures_validated = list(map(lambda r: RibosomeStructure.model_validate(r), structures))
    return { "structures":structures_validated, "count": count }
     
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

class NomenclatureSet(Schema):
    ElongationFactorClass    : list[ElongationFactorClass]
    InitiationFactorClass    : list[InitiationFactorClass]
    CytosolicProteinClass    : list[CytosolicProteinClass ]
    MitochondrialProteinClass: list[MitochondrialProteinClass ]
    CytosolicRNAClass        : list[CytosolicRNAClass ]
    MitochondrialRNAClass    : list[MitochondrialRNAClass ]

@structure_router.get('/list_nomenclature', response=NomenclatureSet)
def polymer_classes_nomenclature(request):
    return {
        "ElongationFactorClass"    : [e.value for e in ElongationFactorClass],
        "InitiationFactorClass"    : [e.value for e in InitiationFactorClass],
        "CytosolicProteinClass"    : [e.value for e in CytosolicProteinClass ],
        "MitochondrialProteinClass": [e.value for e in MitochondrialProteinClass ],
        "CytosolicRNAClass"        : [e.value for e in CytosolicRNAClass ],
        "MitochondrialRNAClass"    : [e.value for e in MitochondrialRNAClass ],
    }


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
