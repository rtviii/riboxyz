from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.ribosome_types.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember

classification_router = Router()
TAG              = "Polymer Classes"

@classification_router.get('/polynucleotide',  tags=[TAG], response=list[RNA])
def polynucleotide_class(request,rna_class:PolynucleotideClass):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA.  """
    agg = []
    rna_class = rna_class.value

    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeAssets(rcsb_id).get_chain_by_polymer_class(rna_class)
        except Exception as e:
            print(e)

        if x is not None:
            agg.append(x.model_dump_json())

    return HttpResponse(agg)

@classification_router.get('/polypeptide',  tags=[TAG], response=list[Protein])
def polypeptide_class(request,protein_class:ProteinClass):
    """All members of the given protein class: cytosolic, mitochondrial proteins """

    agg = []
    protein_class = protein_class.value
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeAssets(rcsb_id).get_chain_by_polymer_class(protein_class)
        except Exception as e:
            print(e)

        if x is not None:
            print("Found a protein class in struct ",protein_class, rcsb_id)
            agg.append(x.model_dump_json())

    return HttpResponse(agg)


@classification_router.get('/lifecyle_factor',  tags=[TAG], response=list[Protein])
def lifecycle_factor_class(request,factor_class:LifecycleFactorClass):
    """All members of the given factor class: initiation, elongation  """
    agg = []
    factor_class = factor_class.value
    for rcsb_id in os.listdir(RIBETL_DATA):
        try:
            x = RibosomeAssets(rcsb_id).get_chain_by_polymer_class(factor_class)
        except Exception as e:
            print(e)
        if x is not None:
            agg.append(x.model_dump_json())

    return HttpResponse(agg)
