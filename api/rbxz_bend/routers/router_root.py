from django.http import HttpResponse
from ninja import NinjaAPI
from ninja.errors import ValidationError
from .router_v0 import v0
from .router_test import test

root = NinjaAPI(docs_url='/docs', openapi_url='/openapi.json',
    title='riboxyz-api',
    description='A programming interface for ribosome data.')


root.add_router('/test', test)
root.add_router('/v0', v0)

# @root.exception_handler(ValidationError)
# def validation_errors(request, exc):
#     return HttpResponse("Invalid input", status=422)
