from django.urls import include, path
from django.conf.urls.static import static
from ribxz_api.settings import STATIC_ROOT, STATIC_URL
from .driver  import root_api
import os

urlpatterns = [
    # path('admin/' , admin.site .urls                   ),
    # path('db/', include('mod_db.urls', 'mod_db')),
    # path('comp/', include('mod_comp.urls', 'mod_comp')),
    path("", root_api.urls),
]+ static(STATIC_URL, document_root=STATIC_ROOT)