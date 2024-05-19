from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl.etl.ribosome_assets import RibosomeAssets
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure

mmcif_router = Router();tag = "pdbx/mmcif_structures"

# E  E - ul4
@mmcif_router.get('/chain',  tags=[tag])
def send_pdbx_mmcif_file(request, rcsb_id:str, auth_asym_id:str):

    rcsb_id = rcsb_id.upper()
    chain_fullpath = os.path.join(os.environ.get('RIBETL_DATA'), rcsb_id, "CHAINS", f'{rcsb_id}_STRAND_{auth_asym_id}.cif')
    
    # Ensure the file exists
    if not os.path.exists(chain_fullpath):
        return HttpResponse(status=404)
    
    # Open and read the file
    with open(chain_fullpath, 'rb') as file:
        file_content = file.read()
    
    # Create the response with the correct MIME type
    response = HttpResponse(file_content, content_type='chemical/x-mmcif')
    response['Content-Disposition'] = f'attachment; filename="{os.path.basename(chain_fullpath)}"'
    
    return response