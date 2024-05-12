from ninja import NinjaAPI, Schema
from ribctl.lib.schema.types_ribosome import CytosolicProteinClass
from routers.router_struct import structure_router
from routers.router_classes import classification_router

# class CytosolicProteinClassSchema(Schema):
#     values: CytosolicProteinClass

# class NinjaWithExtra(NinjaAPI):
#     def get_openapi_schema(self, *args, **kwargs):
#         schema = super().get_openapi_schema(*args, **kwargs)
#         schema["components"]["schemas"]["CytosolicProteinClass"] = CytosolicProteinClassSchema.schema()
#         return schema


root_api = NinjaAPI(docs_url='/', openapi_url='/openapi', title='ribxz API', description='A programming interface for ribosome data.', version='0.1.0' )
root_api.add_router('/structure', structure_router)
root_api.add_router('/polymer_class', classification_router)
