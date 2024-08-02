
from io import BytesIO
import json
import os
from pprint import pprint
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl import RIBETL_DATA
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.libbsite import bsite_ligand
from ribctl.lib.schema.types_binding_site import BindingSite
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, CytosolicProteinClass, PolymerClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure
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
        
