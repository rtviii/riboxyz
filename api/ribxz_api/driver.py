from pprint import pprint
from ninja import NinjaAPI, Schema
from ribctl.lib.schema.types_ribosome import CytosolicProteinClass
from routers.router_struct import structure_router
from routers.router_classes import classification_router
from routers.router_mmcif import mmcif_router


root_api = NinjaAPI(docs_url='/', openapi_url='/openapi', title='ribxz API', description='A programming interface for ribosome data.', version='0.1.0' )
root_api.add_router('/structure', structure_router)
root_api.add_router('/polymer_class', classification_router)
root_api.add_router('/mmcif_structures', mmcif_router)
