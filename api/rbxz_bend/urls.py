from django.urls import include, path
from django.contrib import admin
# TODO: retrofit old urls to the new site.





urlpatterns = [
    # path('admin/' , admin   .site .urls                   ),
    path('db/'    , include('mod_db.urls', 'mod_db')      ),
    path('comp/'  , include('mod_comp.urls', 'mod_comp')  ),
    path('utils/' , include('mod_utils.urls', 'mod_utils')),
    ]
