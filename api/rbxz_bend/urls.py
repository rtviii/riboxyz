from django.urls import include, path
from django.conf.urls.static import static
from rbxz_bend.settings import STATIC_ROOT, STATIC_URL
from .api  import root_api
import os

urlpatterns = [
    # path('admin/' , admin.site .urls                   ),
    # path('db/', include('mod_db.urls', 'mod_db')),
    # path('comp/', include('mod_comp.urls', 'mod_comp')),
    path("", root_api.urls),
]+ static(STATIC_URL, document_root=STATIC_ROOT)