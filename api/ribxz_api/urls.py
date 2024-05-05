from django.urls import include, path
from django.conf.urls.static import static
from ribxz_api.settings import STATIC_ROOT, STATIC_URL
from .driver  import root_api
import os

urlpatterns = [ path("", root_api.urls), ] +  static(STATIC_URL, document_root=STATIC_ROOT)
# urlpatterns = [ root_api.urls ] +  static(STATIC_URL, document_root=STATIC_ROOT)