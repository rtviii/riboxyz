from pprint import pprint
from ninja import NinjaAPI, Schema
from ribctl.lib.schema.types_ribosome import CytosolicProteinClass
from routers.router_struct import structure_router
from routers.router_classes import classification_router
from routers.router_mmcif import mmcif_router


class RibxzAPI(NinjaAPI):
    def get_openapi_schema(self, *a, **kw):
        schema = super().get_openapi_schema(*a, **kw)
        # print("schema", schema)
        struct_list_schema: dict = dict(schema["paths"])["/structures/list"]
        update_ep = {
            "parameters": [
                {
                    "in": "query",
                    "name": "search",
                    "required": False,
                    "schema": {
                        "anyOf": [{"type": "string"}, {"type": "null"}],
                        "title": "Search",
                    },
                },
                {
                    "host_taxa": {
                        "anyOf": [
                            {"items": {"type": "integer"}, "type": "array"},
                            {"type": "null"},
                        ],
                        "title": "Host " "Taxa",
                    },
                    "polymer_classes": {
                        "anyOf": [
                            {
                                "items": {
                                    "anyOf": [
                                        {
                                            "$ref": "#/components/schemas/ribctl__lib__enumunion__UnionEnumMeta__make_union___locals___UnionEnum__1"
                                        },
                                        {
                                            "$ref": "#/components/schemas/ribctl__lib__enumunion__UnionEnumMeta__make_union___locals___UnionEnum__2"
                                        },
                                    ]
                                },
                                "type": "array",
                            },
                            {"type": "null"},
                        ],
                        "title": "Polymer " "Classes",
                    },
                    "resolution": {
                        "anyOf": [
                            {
                                "maxItems": 2,
                                "minItems": 2,
                                "prefixItems": [
                                    {"anyOf": [{"type": "number"}, {"type": "null"}]},
                                    {"anyOf": [{"type": "number"}, {"type": "null"}]},
                                ],
                                "type": "array",
                            },
                            {"type": "null"},
                        ],
                        "title": "Resolution",
                    },
                    "source_taxa": {
                        "anyOf": [
                            {"items": {"type": "integer"}, "type": "array"},
                            {"type": "null"},
                        ],
                        "title": "Source " "Taxa",
                    },
                    "year": {
                        "anyOf": [
                            {
                                "maxItems": 2,
                                "minItems": 2,
                                "prefixItems": [
                                    {"anyOf": [{"type": "integer"}, {"type": "null"}]},
                                    {"anyOf": [{"type": "integer"}, {"type": "null"}]},
                                ],
                                "type": "array",
                            },
                            {"type": "null"},
                        ],
                        "title": "Year",
                    },
                },
            ]
        }
        struct_list_schema.update(update_ep)
        schema.update(struct_list_schema)

        return schema


api = RibxzAPI()

root_api = RibxzAPI(
    docs_url="/",
    openapi_url="/openapi",
    title="ribxz API",
    description="A programming interface for ribosome data.",
    version="0.1.0",
)
root_api.add_router("/structures", structure_router)
root_api.add_router("/polymers", classification_router)
root_api.add_router("/mmcif", mmcif_router)
