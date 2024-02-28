from ninja import NinjaAPI
from .api_v0_router import apiv0
from .api_structure_router import structure_router

root_api = NinjaAPI(docs_url='/docs', openapi_url='/openapi.json', title='riboxyz-api', description='A programming interface for ribosome data.')
root_api.add_router('/', apiv0)
root_api.add_router('/structure', structure_router)