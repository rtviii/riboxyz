from django.conf import settings
from django.urls import include, path
from django.contrib import admin
from django.conf.urls.static import static
# TODO: retrofit old urls to the new site.





urlpatterns = [
    # path('admin/' , admin   .site .urls                   ),
    path('db/'    , include('mod_db.urls', 'mod_db')      ),
    path('comp/'  , include('mod_comp.urls', 'mod_comp')  ),
    path('utils/' , include('mod_utils.urls', 'mod_utils')),
    ]

urlpatterns += static(settings.STATIC_URL, document_root=settings.STATIC_ROOT)