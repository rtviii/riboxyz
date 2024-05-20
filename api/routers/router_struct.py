import json
from pprint import pprint
import typing
from django.http import  JsonResponse, HttpResponseServerError
from ninja import Router, Schema
from api.ribxz_api.db_queries import dbqueries
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import  CytosolicProteinClass, CytosolicRNAClass, ElongationFactorClass, InitiationFactorClass, LifecycleFactorClass, MitochondrialProteinClass, MitochondrialRNAClass, PolymerClass, PolynucleotideClass, PolypeptideClass, ProteinClass, RibosomeStructure, tRNA
from ribctl.lib.libtax import Taxid, ncbi

structure_router = Router()
TAG              = "Structure"

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
        pprint(ptc)
    except Exception as e:
        return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))
    # except Exception as e:
    #     return HttpResponseServerError("Failed to find structure profile {}:\n\n{}".format(rcsb_id, e))

@structure_router.get('/list_structures', response=list[RibosomeStructure], tags=[TAG])
def list_structures(request):
    structs_response = dbqueries.list_structs()
    structures       = list(map(lambda r: RibosomeStructure.model_validate(r), structs_response))
    return structures



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
