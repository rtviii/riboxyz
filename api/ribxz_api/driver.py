from pprint import pprint
from ninja import NinjaAPI, Schema
from routers.router_loci import router_loci
from routers.router_struct import router_structures
from routers.router_polymers import router_polymers
from routers.router_mmcif import router_mmcif
from routers.router_lig import router_lig



root_api = NinjaAPI(
    docs_url    = "/",
    openapi_url = "/openapi",
    title       = "ribxz API",
    description = "A programming interface for ribosome data.",
    version     = "0.1.0",
)
root_api.add_router("/structures", router_structures)
root_api.add_router("/polymers", router_polymers)
root_api.add_router("/ligands", router_lig)
root_api.add_router("/loci", router_loci)
root_api.add_router("/mmcif", router_mmcif)
