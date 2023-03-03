from django.conf import settings
from django.urls import include, path
from django.contrib import admin
from django.conf.urls.static import static
from rbxz_bend.api import router
# TODO: retrofit old urls to the new site.
from ninja import NinjaAPI
api = NinjaAPI()
api.add_router('/blog/',router)




urlpatterns = [
    # path('admin/' , admin   .site .urls                   ),
    path('db/'    , include('mod_db.urls', 'mod_db')      ),
    path('comp/'  , include('mod_comp.urls', 'mod_comp')  ),
    path('utils/' , include('mod_utils.urls', 'mod_utils')),
    path('api/',api.urls)
    ]

# urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)