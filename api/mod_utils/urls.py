from django.urls import path
from .views import * 

#TODO:
# - map old endpoints to new ones
# - serve pfam mappings to protein families
# - "last updated" endpoint for fend serving : date, some new structure ids.

urlpatterns = [
    path('hello/' , hello    ),
    path('struct_assets/' , struct_assets    ),
    path('pull_struct_pdb/', pull_struct_pdb),
]
app_name = 'mod_utils'