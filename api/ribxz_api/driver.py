from pprint import pprint
from ninja import NinjaAPI, Schema
from routers.router_struct import structure_router
from routers.router_classes import classification_router
from routers.router_mmcif import mmcif_router
from routers.router_lig import router_lig



root_api = NinjaAPI(
    docs_url    = "/",
    openapi_url = "/openapi",
    title       = "ribxz API",
    description = "A programming interface for ribosome data.",
    version     = "0.1.0",
)
root_api.add_router("/structures", structure_router)
root_api.add_router("/polymers", classification_router)
root_api.add_router("/mmcif", mmcif_router)
root_api.add_router("/ligand", router_lig)
