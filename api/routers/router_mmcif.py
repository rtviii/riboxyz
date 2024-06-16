from io import BytesIO
import json
import os
from django.http import HttpResponse, JsonResponse, HttpResponseServerError
from ninja import Router
from ribctl.etl.etl_assets_ops import RibosomeOps, Structure
from ribctl.lib.schema.types_ribosome import RNA, LifecycleFactorClass, MitochondrialProteinClass, Polymer, PolymerClass, CytosolicProteinClass, PolynucleotideClass, PolypeptideClass, Protein, ProteinClass, RibosomeStructure

mmcif_router = Router();
tag = "mmcif"

# E  E - ul4
@mmcif_router.get('/polymer',  tags=[tag])
def polymer(request, rcsb_id:str, auth_asym_id:str):
    rcsb_id        = rcsb_id.upper()
    RO             = RibosomeOps(rcsb_id)
    chain_fullpath = os.path.join(RO.paths.chains_dir, f'{rcsb_id}_{auth_asym_id}.cif')
    
    print( chain_fullpath)
    # Ensure the file exists
    if not os.path.exists(chain_fullpath):
        return HttpResponse(status=404)
    
    # Open and read the file
    with open(chain_fullpath, 'rb') as file:
        file_content = file.read()
    print("GOT FILE AT APTH" , chain_fullpath)
    
    # Create the response with the correct MIME type
    response = HttpResponse(file_content, content_type='chemical/x-mmcif')
    response['Content-Disposition'] = f'attachment; filename="{os.path.basename(chain_fullpath)}"'
    
    return response