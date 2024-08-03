
from io import BytesIO
import json
import os
from pprint import pprint
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.libbsite import bsite_ligand, init_transpose_ligand
from ribctl.lib.schema.types_binding_site import BindingSite, LigandTransposition
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolynucleotideClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember

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
def lig_transpose(request, source_structure:str, target_structure:str, chemical_id:str, radius:int=5):
    """All members of the given RNA class: small and large subunit, cytosolic and mitochondrial RNA; tRNA.  """
    bsite_path      = RibosomeOps(source_structure).paths.binding_site(chemical_id)
    prediction_path = RibosomeOps(target_structure).paths.binding_site_prediction(chemical_id, source_structure)

    if not os.path.exists(bsite_path):
        try:
            bsite = bsite_ligand(chemical_id, source_structure, radius, save=True)
            if not os.path.exists(prediction_path):
                prediction = init_transpose_ligand(RibosomeOps(target_structure).profile(),bsite)
                with open(prediction_path, 'w') as f:
                    json.dump(prediction.model_dump(), f)
                    print("Saved {}".format(prediction_path))
                return JsonResponse(prediction.model_dump(), safe=False)
        except Exception as e:
            return HttpResponseServerError(e)
    else:
        with open(bsite_path, 'r') as f:
            data = json.load(f)
            bsite = BindingSite.model_validate(data)
        prediction = init_transpose_ligand(RibosomeOps(target_structure).profile(),bsite)

        with open(prediction_path, 'w') as f:
            json.dump(prediction.model_dump(), f)
            print("Saved {}".format(prediction_path))
        return JsonResponse(prediction.model_dump(), safe=False)

