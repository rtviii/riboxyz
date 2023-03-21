from ninja import NinjaAPI
from .router_v0 import v0
from .router_test import test

root = NinjaAPI(docs_url='/docs', openapi_url='/openapi.json',
    title='riboxyz-api',
    description='A programming interface for ribosome data.')

root.add_router('/test', test)
root.add_router('/v0', v0)
