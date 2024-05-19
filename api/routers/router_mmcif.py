from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure
from schema.v0 import BanClassMetadata, ExogenousRNAByStruct,LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember

mmcif_router = Router();tag = "pdbx/mmcif_structures"

@mmcif_router.get('/mmcif_chain',  tags=[tag])
def send_pdbx_mmcif_file(request, rcsb_id:str, auth_asym_id:str):
    print("Got ", rcsb_id, auth_asym_id)
    print(request)
    
    # full_path = os.path.join(settings.MEDIA_ROOT, file_path)
    
    # Ensure the file exists
    if not os.path.exists(full_path):
        return HttpResponse(status=404)
    
    # Open and read the file
    with open(full_path, 'rb') as file:
        file_content = file.read()
    
    # Create the response with the correct MIME type
    response = HttpResponse(file_content, content_type='chemical/x-mmcif')
    response['Content-Disposition'] = f'attachment; filename="{os.path.basename("3j7z.cif")}"'
    
    return response