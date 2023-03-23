from django.http import HttpResponse
from ninja import NinjaAPI
from .router_v0 import v0
from .router_test import test
from ninja.errors import ValidationError

root = NinjaAPI(docs_url='/', openapi_url='/openapi.json',
    title='riboxyz-api',
    description='A programming interface for ribosome data.')

# @root.exception_handler(ValidationError)
# def validation_errors(request, exc):
#     return HttpResponse("some validation failed", status=200)

@root.exception_handler(ValidationError)
def custom_validation_errors(request, exc):
    print(exc.errors)  # <--------------------- !!!!
    return root.create_response(request, {"detail": exc.errors}, status=422)

root.add_router('/test', test)
root.add_router('/v0', v0)

