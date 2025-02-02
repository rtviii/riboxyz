import json
import os
from pprint import pprint
from django.http import JsonResponse
from ninja import Router
from neo4j_ribosome.db_lib_reader import PolymersFilterParams
from ribctl import RIBETL_DATA
from ribctl.ribosome_ops import RibosomeOps
from ribctl.lib.schema.types_ribosome import RNA, Polymer, Protein
from ribctl.lib.types.polymer import (
    LifecycleFactorClass,
    PolynucleotideClass,
    PolynucleotideClass,
    PolypeptideClass,
)

from neo4j_ribosome.db_lib_reader import (
    PolymersFilterParams,
    StructureFilterParams,
    dbqueries,
)

router_polymers = Router()
TAG = "Polymer Classes"


@router_polymers.get("/polynucleotide", tags=[TAG], response=list[RNA])
def polynucleotide_class(request, rna_class: PolynucleotideClass):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA."""
    agg = []
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeOps(rcsb_id).get_poly_by_polyclass(rna_class)
        except Exception as e:
            print(e)

        if x is not None:
            agg.append(x.model_dump())

    return JsonResponse(agg, safe=False)


@router_polymers.get("/polypeptide", tags=[TAG], response=list[Protein])
def polypeptide_class(request, protein_class: PolypeptideClass):
    """All members of the given protein class: cytosolic, mitochondrial proteins"""

    agg = []
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeOps(rcsb_id).get_poly_by_polyclass(protein_class)
        except Exception as e:
            print(e)

        if x is not None:
            agg.append(x.model_dump())

    return JsonResponse(agg, safe=False)


@router_polymers.get("/lifecyle_factor", tags=[TAG], response=list[Protein])
def lifecycle_factor_class(request, factor_class: LifecycleFactorClass):
    """All members of the given factor class: initiation, elongation"""
    agg = []
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeOps(rcsb_id).get_poly_by_polyclass(factor_class)
        except Exception as e:
            print(e)
        if x is not None:
            agg.append(x.model_dump())

    return JsonResponse(agg, safe=False)


@router_polymers.post("/list_polymers", response=dict, tags=[TAG])
def list_polymers(request, filters: PolymersFilterParams):
    parsed_filters = PolymersFilterParams(**json.loads(request.body))
    polymers, next_cursor, total_polymers_count, total_structures_count = ( dbqueries.list_polymers_filtered(parsed_filters) )
    polymers_validated = [Polymer.model_validate(p) for p in polymers]
    return {
        "polymers": polymers_validated,
        "next_cursor": next_cursor,
        "total_polymers_count": total_polymers_count,
        "total_structures_count": total_structures_count,
    }