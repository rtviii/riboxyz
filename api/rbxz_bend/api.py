from ninja import NinjaAPI

root_api = NinjaAPI(docs_url='/docs', openapi_url='/openapi.json', title='riboxyz-api', description='A programming interface for ribosome data.')

@root_api.get("/hello")
def hello(request):
    return "Hello world"