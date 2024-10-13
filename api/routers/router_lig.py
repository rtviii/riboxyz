
from io import BytesIO
import json
import os
from pprint import pprint
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.libbsite import bsite_ligand, bsite_transpose, lig_get_chemical_categories
from ribctl.lib.schema.types_binding_site import BindingSite, LigandTransposition
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolynucleotideClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructureMetadata

router_lig = Router()
TAG                   = "Ligands, Antibitics & Small Molecules"

@router_lig.get('/binding_pocket',  tags=[TAG], response=BindingSite)
def lig_nbhd(request, source_structure:str, chemical_id:str, radius:int=5):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA.  """
    path = RibosomeOps(source_structure).paths.binding_site(chemical_id)
    if not os.path.exists(path):
        try:
            bsite = bsite_ligand(chemical_id, source_structure, radius)
            pprint("produced")
            pprint(bsite)
            return JsonResponse(bsite.model_dump(), safe=False)
        except Exception as e:
            return HttpResponseServerError(e)
       
    with open(path, 'r') as f:
        return JsonResponse(json.load(f), safe=False)
        
@router_lig.get('/transpose',  tags=[TAG], response=LigandTransposition)
def lig_transpose(request, source_structure:str, target_structure:str, chemical_id:str, radius:float):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA.  """
    bsite_path      = RibosomeOps(source_structure).paths.binding_site(chemical_id)
    prediction_path = RibosomeOps(target_structure).paths.binding_site_prediction(chemical_id, source_structure)

    bsite           = bsite_ligand(chemical_id, source_structure, radius)
    prediction      = bsite_transpose(source_structure,target_structure,bsite)

    return JsonResponse(prediction.model_dump(), safe=False)

@router_lig.get('/chemical_classification',  tags=[TAG], response=dict)
def lig_chemical_categories(request,):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA.  """
    return JsonResponse(lig_get_chemical_categories(), safe=False)



