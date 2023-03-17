import sys
import os
from django.conf import settings
from django.urls import include, path
from django.contrib import admin
from django.conf.urls.static import static
from ninja import NinjaAPI
from rbxz_bend.settings import STATIC_ROOT, STATIC_URL
import asyncio
import random
import typing
from venv import logger
from ninja import Router
from rbxz_bend.settings import get_logger
from ribctl.lib.struct_rcsb_api import current_rcsb_structs
from ribctl.lib.types.types_ribosome_assets import RibosomeAssets
from ribctl.db.ribosomexyz import riboxyzDB
from ribctl.lib.types.types_polymer import RNAClass
from ribctl.lib.types.types_ribosome import ExogenousRNAByStruct, ProteinClass, RibosomeStructure
from schema.v0 import BanClassMetadata, LigandInstance, LigandlikeInstance, NeoStruct, NomenclatureClass, NomenclatureClassMember
from .routers import router_root
import os

urlpatterns = [
    # path('admin/' , admin.site .urls                   ),
    path('db/', include('mod_db.urls', 'mod_db')),
    path('comp/', include('mod_comp.urls', 'mod_comp')),
    path('', router_root.root.urls),
]+ static(STATIC_URL, document_root=STATIC_ROOT)