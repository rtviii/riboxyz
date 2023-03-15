import sys
import os
from django.conf import settings
from django.urls import include, path
from django.contrib import admin
from django.conf.urls.static import static
from ninja import NinjaAPI
from .api import ninja_api_router

api = NinjaAPI(
    title='riboxyz-API',
)
api.add_router('', ninja_api_router)


urlpatterns = [
    # path('admin/' , admin   .site .urls                   ),
    path('db/', include('mod_db.urls', 'mod_db')),
    path('comp/', include('mod_comp.urls', 'mod_comp')),
    path('utils/', include('mod_utils.urls', 'mod_utils')),
    path('api/', api.urls)
]

# urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)

